/*----------------------------------------------------------------------
| This file contains routines and global variables that are common for
| all operating systems the program has been ported to.  It is included
| in one of the source code files of each port.  See common.h for the
| common #defines and common routine definitions.
|
| Commona contains information used only during setup
| Commonb contains information used only during execution
| Commonc contains information used during setup and execution
|
| Copyright 1995-2019 Mersenne Research, Inc.  All rights reserved
+---------------------------------------------------------------------*/

/* Routine to eliminate odd puctuation characters from user ID */
/* and computer ID */

void sanitizeString (
        char    *p)
{
        int     i;
        for (i = (int) strlen (p); i > 0 && isspace (p[i-1]); i--) p[i-1] = 0;
        while (*p) {
                if (!IsCharAlphaNumeric (*p) &&
                    *p != '.' && *p != '-' && *p != '_')
                        *p = '_';
                p++;
        }
}

/* Create a status report message from the work-to-do file */

#define STAT0 "Below is a report on the work you have queued and any expected completion dates.\n"
#define STAT1 "The chance that one of the %d exponents you are testing will yield a %sprime is about 1 in %lld. "
#define STAT1a "The chance that the exponent you are testing will yield a %sprime is about 1 in %lld. "
#define STAT3 "No work queued up.\n"

void rangeStatusMessage (
        char    *buf,
        unsigned int buflen)            /* Originally coded for a 2000 character buffer */
{
        unsigned int tnum, ll_and_prp_cnt, lines_per_worker;
        int     mersennes;              /* TRUE if only testing Mersenne numbers */
        double  prob, est;
        char    *orig_buf;

/* Just in case the user hand added work to the worktodo file, reread it */
/* now if the worker threads and communication threads are not active. */

        if (! WORKER_THREADS_ACTIVE && !COMMUNICATION_THREAD) readIniFiles ();

/* Init.  Default is 32 lines in a 2000 character buffer */

        lines_per_worker = (unsigned int) IniGetInt (INI_FILE, "StatusLines", buflen / 62) / NUM_WORKER_THREADS;
        if (lines_per_worker < 3) lines_per_worker = 3;
        orig_buf = buf;
        ll_and_prp_cnt = 0;
        prob = 0.0;
        mersennes = TRUE;
        strcpy (buf, STAT0);
        buf += strlen (buf);

/* Loop over all worker threads */

        for (tnum = 0; tnum < NUM_WORKER_THREADS; tnum++) {
            struct work_unit *w;
            unsigned int lines_output;
            int truncated_status_msg;

/* Init line formatting info */

            lines_output = 0;
            truncated_status_msg = FALSE;

/* Output thread id */

            if (NUM_WORKER_THREADS > 1) {
                sprintf (buf, "[Worker thread #%d]\n", tnum+1);
                buf += strlen (buf);
                lines_output++;
            }

/* Loop over all work units */

            w = NULL;
            est = 0.0;
            for ( ; ; ) {
                time_t  this_time;
                char    timebuf[80];
                unsigned int bits;

/* Read the next line of the work file */

                w = getNextWorkToDoLine (tnum, w, SHORT_TERM_USE);
                if (w == NULL) break;
                if (w->work_type == WORK_NONE) continue;

/* Keep track of whether we are only testing Mersenne numbers */

                if (w->k != 1.0 || w->b != 2 || w->c != -1 || w->known_factors != NULL) mersennes = FALSE;

/* If testing then adjust our probabilities */
/* This assumes our error rate is roughly 1.8% */

                bits = (unsigned int) w->sieve_depth;
                if (bits < 32) bits = 32;
                if (w->work_type == WORK_TEST) {
                        ll_and_prp_cnt++;
                        prob += (bits - 1) * 1.733 * (w->pminus1ed ? 1.04 : 1.0) / (_log2(w->k) + _log2(w->b) * w->n);
                }
                if (w->work_type == WORK_DBLCHK) {
                        ll_and_prp_cnt++;
                        prob += (bits - 1) * 1.733 * ERROR_RATE * (w->pminus1ed ? 1.04 : 1.0) / (_log2(w->k) + _log2(w->b) * w->n);
                }
                if (w->work_type == WORK_PRP) {
                        ll_and_prp_cnt++;
                        if (!w->prp_dblchk)
                                prob += (bits - 1) * 1.733 * (w->pminus1ed ? 1.04 : 1.0) / (_log2(w->k) + _log2(w->b) * w->n);
                        else
                                prob += (bits - 1) * 1.733 * PRP_ERROR_RATE * (w->pminus1ed ? 1.04 : 1.0) / (_log2(w->k) + _log2(w->b) * w->n);
                }

/* Adjust our time estimate */

                est += work_estimate (tnum, w);

/* Stop adding worktodo lines if buffer is full.  We must still loop */
/* through the worktodo lines to decrement the in-use counters. */

                if ((unsigned int) (buf - orig_buf) >= buflen - 200 ||
                    lines_output >= lines_per_worker-1) {
                        if (! truncated_status_msg) {
                                strcpy (buf, "More...\n");
                                buf += strlen (buf);
                                truncated_status_msg = TRUE;
                        }
                        continue;
                }

/* Add the exponent to the output message */

                gw_as_string (buf, w->k, w->b, w->n, w->c);
                buf += strlen (buf);
                if (w->work_type == WORK_PRP && w->known_factors) {
                        strcpy (buf, "/known_factors");
                        buf += strlen (buf);
                }
                strcpy (buf, ", ");
                buf += strlen (buf);

                if (w->work_type == WORK_ECM)
                        sprintf (buf, "ECM %d curve%s B1=%.0f",
                                 w->curves_to_do,
                                 w->curves_to_do == 1 ? "" : "s",
                                 w->B1);
                else if (w->work_type == WORK_PMINUS1)
                        sprintf (buf, "P-1 B1=%.0f", w->B1);
                else if (w->work_type == WORK_FACTOR)
                        sprintf (buf, "factor from 2^%d to 2^%d",
                                 (int) w->sieve_depth, (int) w->factor_to);
                else
                        strcpy (buf, w->work_type == WORK_PFACTOR ? "P-1" :
                                     w->work_type == WORK_TEST ||
                                     w->work_type == WORK_ADVANCEDTEST ? "Lucas-Lehmer test" :
                                     w->work_type == WORK_DBLCHK ? "Double-check" :
                                     /* w->work_type == WORK_PRP */ "PRP");
                buf += strlen (buf);

                time (&this_time);
                if (est + (double) this_time < 2147483640.0) {
                        this_time += (long) est;
                        strcpy (timebuf, ctime (&this_time));
                        safe_strcpy (timebuf+16, timebuf+19);
                } else
                        strcpy (timebuf, "after Jan 19 2038\n");
                sprintf (buf, ", %s", timebuf);
                buf += strlen (buf);
                lines_output++;
            }

/* Format more of the message */

            if (est == 0.0 && ! truncated_status_msg) {
                strcpy (buf, STAT3);
                buf += strlen (buf);
            }
        }

/* Print message estimating our probability of success */

        if (ll_and_prp_cnt == 1)
                sprintf (buf+strlen(buf), STAT1a, mersennes ? "Mersenne " : "", (long long) (1.0 / prob));
        if (ll_and_prp_cnt > 1)
                sprintf (buf+strlen(buf), STAT1, ll_and_prp_cnt, mersennes ? "Mersenne " : "", (long long) (1.0 / prob));
}



/* Return the suggested minimum number of cores that should be used for a work preference. */
/* Used in the Worker Windows dialog box. */

int min_cores_for_work_pref (
        int     work_pref)
{
        int     cores;

// Default minimum number of cores is 1.

        cores = 1;

// If LL or PRP testing 100M digit numbers, use at least 4 cores (or all cores)

        if (work_pref == PRIMENET_WP_LL_100M || work_pref == PRIMENET_WP_PRP_100M) {
                if (NUM_CPUS < 8) cores = NUM_CPUS;
                else cores = 4;
        }

// If we aren't using the computer 24 hours a day, then scale the minimum number of cores up

        cores = cores * 24 / CPU_HOURS;
        if (cores > (int) NUM_CPUS) cores = NUM_CPUS;

// Return the minimum number of cores

        return (cores);
}


/**********CLEANUP******************/
/* TODO make opendir & readdir generics for mac/windows */

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
/**********CLEANUP******************/


/* Create a status report message about Backup/Restore files in working directory */

#define BACKUP_NONE         "No Backup/Restore files (*.bu) were found in %s.\n"
#define BACKUP_CWD_ERROR    "Unable to read working directory (%s).\n"
#define BACKUP_PARSE_ERROR  "Unable to parse (%s).\n"

void restoreStatusMessage (
        char    *buf,
        unsigned int buflen)            /* Originally coded for a 1000 character buffer */
{
        unsigned int ll_and_prp_cnt, ecm_and_pm1_cnt;
        char    *orig_buf;

/* TODO augment by checking how many of these are mentioned in worktodo. */

/* Init.  Default is 16 lines in a 1000 character buffer */
        orig_buf = buf;

/* TODO add ll_and_prp_cnt and ecm_and_pm1_cnt */
/* TODO get wDir name so it can be added to status message */

// #include <unistd.h>
//    if (getcwd(cwd, sizeof(cwd)) != NULL) {
//      printf("Current working dir: %s\n", cwd);

        struct dirent *dp;
        DIR *dfd;
        if ((dfd = opendir(".") == NULL)
        {
                sprintf(buf, BACKUP_CWD_ERROR, "TODO");
                buf += strlen(buf);
                return (FALSE);
        }

        char filename_qfd[100] ;
        while ((dp = readdir(dfd)) != NULL)
        {
                sprintf( filename_qfd , "%s",dir,dp->d_name) ;
                if( dp->d_type != DT_REG )
                {
                        continue ;
                }

                char *ext = strrchr(dp->d_name, '.');
                if (ext && strcmp(ext, ".bu")) {
                        // Load File Status
                        struct work_unit *w;
                        if (!restoreStatus(dp->d_name, work_unit))
                        {
                                snprintf(buf, buflen, BACKUP_PARSE_ERROR, dp->d_name);
                                buflen -= strlen(buf);
                                buf += strlen(buf);
                        }
                }
        }
}

int restoreFileStatusMessage (
        char    *filename,
        struct work_unit *w)
{
        int        fd;
        unsigned long file_magicnum;
        unsigned long version;
        double  k;
        unsigned long b, n;
        long    c;
        char    pad;
        char    stage[11];
        double  pct_complete;
        unsigned long tmp;

// See pct_complete_from_savefile
// readLLSaveFile
// read_header
// ecm_restore
// pm1_restore

        fd = _open (filename, _O_BINARY | _O_RDONLY);
        if (fd <= 0) return (FALSE);

/* Load the file magicnum */

        // read_magicnum & read_header don't return the read values.
        // reusing portions of that code here.

        _lseek(fd, 0, SEEK_SET);
        if (!read_long (fd, &file_magicnum, NULL)) return (FALSE);

/* Load the rest of the common file header */
        if (!read_long (fd, version, NULL)) return (FALSE);

        if (!read_double (fd, &k, NULL)) return (FALSE);
        if (!read_long (fd, &b, NULL)) return (FALSE);
        if (!read_long (fd, &n, NULL)) return (FALSE);
        if (!read_slong (fd, &c, NULL)) return (FALSE);

/* Update workunit fields */

        w->k = k;
        w->b = b;
        w->n = n;
        w->c = c;

/* Can read_header to validate and set other fields */

        if (!read_header(fd, &trash, w, NULL)) return (FALSE);

        switch (file_magicnum) {
        case: ECM_MAGICNUM
                w->work_type =  WORK_ECM;
                if (version != ECM_VERSION) return (FALSE);
                if (! read_long (fd, &tmp, NULL)) goto readerr;
                // TODO Where to write stage???
                //*stage = (int) tmp; 
                if (! read_long (fd, w->curve, &sum)) goto readerr;
                if (! read_long (fd, w->curve, &sum)) goto readerr;


                break;
        case: PM1_MAGICNUM
                w->work_type =  WORK_PMINUS1;
                if (version != PM1_VERSION) return (FALSE);
                break;
        case: LL_MAGICNUM
                w->work_type =  WORK_TEST;
                if (version != LL_VERSION) return (FALSE);
                break;
        case: PRP_MAGICNUM
                w->work_type =  WORK_PRP;
                if (version != PRP_VERSION) return (FALSE);
                break;
        case: FACTOR_MAGICNUM
                w->work_type =  WORK_FACTOR;
                if(version != FACTOR_MAGICNUM) return (FALSE);
                break;
        default:
            return (FALSE);
        }

/* Return success */
        return (TRUE);

}



        if (!read_magicnum (fd, PRP_MAGICNUM)) goto err;
        if (!read_header (fd, &version, w, &filesum)) goto err;
        if (version == 0 || version > PRP_VERSION) goto err;



        if (! read_magicnum (fd, PM1_MAGICNUM)) goto readerr;
        if (! read_header (fd, &version, w, &filesum)) goto readerr;
        if (version != 1 && version != 2) goto readerr;

/* Read the file data */

        if (! read_long (fd, &pm1data->stage, &sum)) goto readerr;
        if (! read_longlong (fd, &pm1data->B_done, &sum)) goto readerr;
        if (! read_longlong (fd, &pm1data->B, &sum)) goto readerr;
        if (! read_longlong (fd, &pm1data->C_done, &sum)) goto readerr;
        if (! read_longlong (fd, &pm1data->C_start, &sum)) goto readerr;
        if (! read_longlong (fd, &pm1data->C, &sum)) goto readerr;
        if (! read_longlong (fd, processed, &sum)) goto readerr;
        if (! read_long (fd, &pm1data->D, &sum)) goto readerr;
        if (! read_long (fd, &pm1data->E, &sum)) goto readerr;
        if (! read_long (fd, &pm1data->rels_done, &sum)) goto readerr;
        if (! read_long (fd, &pm1data->bitarray_len, &sum)) goto readerr;


        if ((unsigned int) (buf - orig_buf) >= buflen - 200 ||
            lines_output >= lines_per_worker-1) {
                if (! truncated_status_msg) {
                        strcpy (buf, "More...\n");
                        buf += strlen (buf);
                        truncated_status_msg = TRUE;
                }
                continue;
        }

        if (ll_and_prp_cnt == 1)
                sprintf (buf+strlen(buf), STAT1a, mersennes ? "Mersenne " : "", (long long) (1.0 / prob));
        if (ll_and_prp_cnt > 1)
                sprintf (buf+strlen(buf), STAT1, ll_and_prp_cnt, mersennes ? "Mersenne " : "", (long long) (1.0 / prob));
}





/* Return the suggested minimum number of cores that should be used for a work preference. */
/* Used in the Worker Windows dialog box. */

int min_cores_for_work_pref (
        int     work_pref)
{
        int     cores;

// Default minimum number of cores is 1.

        cores = 1;

// If LL or PRP testing 100M digit numbers, use at least 4 cores (or all cores)

        if (work_pref == PRIMENET_WP_LL_100M || work_pref == PRIMENET_WP_PRP_100M) {
                if (NUM_CPUS < 8) cores = NUM_CPUS;
                else cores = 4;
        }

// If we aren't using the computer 24 hours a day, then scale the minimum number of cores up

        cores = cores * 24 / CPU_HOURS;
        if (cores > (int) NUM_CPUS) cores = NUM_CPUS;

// Return the minimum number of cores

        return (cores);
}

