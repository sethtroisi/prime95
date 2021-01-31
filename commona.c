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


/**********CLEANUP******************/
/* TODO make opendir, readdir, getcwd generics for mac/windows */

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
/**********CLEANUP******************/


/* Create a status report message about backup/restore files in working directory */

#define MAX_BACKUP_FILES    100
#define BACKUP_CWD          "Status of files in '%s'.\n"
#define BACKUP_CWD_ERROR    "Unable to read working directory.\n"
#define BACKUP_STATUS       "Backup %-16s | %s.\n"
#define BACKUP_NONE         "No Backup files (*.bu) were found in %s.\n"
#define BACKUP_PARSE_ERROR  "Unable to parse (%s).\n"

<<<<<<< Updated upstream
int restoreWorkUnitFromFile (
        char    *filename,
        struct work_unit *w,
        pm1handle *pm1)
=======
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
        if ((dfd = opendir(".")) == NULL)
        {
                sprintf(buf, BACKUP_CWD_ERROR, "TODO");
                buf += strlen(buf);
                return;
        }

        // TODO what is qfd
        //char filename_qfd[100] ;
        while ((dp = readdir(dfd)) != NULL)
        {
                //sprintf( filename_qfd , "%s/%s", dir , dp->d_name) ;
                if( dp->d_type != DT_REG )
                {
                        continue ;
                }

                char *ext = strrchr(dp->d_name, '.');
                if (ext && strcmp(ext, ".bu")) {
                        // Load File into work_unit.
                        struct work_unit w;
                        pm1handle pm1;
                        if (!restoreWorkUnitFromFile(dp->d_name, &w, &pm1))
                        {
                                snprintf(buf, buflen, BACKUP_PARSE_ERROR, dp->d_name);
                                buflen -= strlen(buf);
                                buf += strlen(buf);
                        }
                }
        }
}

int restoreWorkUnitFromFile (
        char    *filename,
        struct work_unit *w,
        struct pm1handle *pm1)
>>>>>>> Stashed changes
{
        int        fd;
        unsigned long file_magicnum;
        unsigned long version;
<<<<<<< Updated upstream
=======
        double  k;
        unsigned long b, n;
        signed long    c;
>>>>>>> Stashed changes
        char    pad;
        char    stage[11];
        double  pct_complete;
        unsigned long tmp;

        w->work_type = WORK_NONE;

        fd = _open (filename, _O_BINARY | _O_RDONLY);
        if (fd <= 0) goto readerr;

/* Load the file magicnum. */

        // read_magicnum & read_header don't return the read values.
        // reusing portions of that code here.

        _lseek(fd, 0, SEEK_SET);
        if (!read_long (fd, &file_magicnum, NULL)) goto readerr;

/* Load the rest of the common file header into the workunit. */

        if (!read_long (fd, &version, NULL)) goto readerr;

        if (!read_double (fd, &w->k, NULL)) goto readerr;
        if (!read_long (fd, &w->b, NULL)) goto readerr;
        if (!read_long (fd, &w->n, NULL)) goto readerr;
        if (!read_slong (fd, &w->c, NULL)) goto readerr;

/* Call read_header to validate and set some other fields */

        if (!read_header(fd, &tmp, w, NULL)) goto readerr;

/* Load work type specific data. */

// TODO move _MAGICNUM and _VERSION to header.
        switch (file_magicnum) {
<<<<<<< Updated upstream
        case 0x1725bcd9: // ECM_MAGICNUM:
                //if (version != ECM_VERSION) goto readerr;
                if (version != 1) goto readerr;

                // duplicated from ecm.c ecm_restore (minus reading data)
                w->work_type =  WORK_ECM;

                if (! read_long (fd, &pm1->stage, NULL)) goto readerr;
                if (! read_long (fd, &tmp, NULL)) goto readerr;
                w->curves_to_do = (int) tmp;
                // Sigma
                if (! read_double (fd, &w->curve, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->B, NULL)) goto readerr;
                // Stage 1 current P
                if (! read_longlong (fd, &pm1->B_done, NULL)) goto readerr;
                // Stage 2 current P
                if (! read_longlong (fd, &pm1->C_done, NULL)) goto readerr;
=======
        case: ECM_MAGICNUM
                // duplicated from ecm.c (minus reading data)
                w->work_type =  WORK_ECM;
                if (version != ECM_VERSION) return (FALSE);

                if (! read_long (fd, &tmp, NULL)) return (FALSE)
                w->stage[0] = (char) tmp;
                if (! read_long (fd, w->curves_to_do, NULL)) goto readerr;
                if (! read_double (fd, w->curve, NULL)) goto readerr;
                if (! read_longlong (fd, w-> B, NULL)) goto readerr;
                if (! read_longlong (fd, w-> B_PROCESSED, NULL)) goto readerr;
                if (! read_longlong (fd, w-> C_PROCESSED, NULL)) goto readerr;
>>>>>>> Stashed changes
                break;
        case 0x317a394b: // PM1_MAGICNUM:
                //if (version != PM1_VERSION) goto readerr;
                if (version != 2) goto readerr;

                // duplicated from ecm.c pm1_restore (minus reading data)
                w->work_type =  WORK_PMINUS1;
<<<<<<< Updated upstream

                if (! read_long (fd, &pm1->stage, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->B_done, NULL)) goto readerr;
=======
                if (version != PM1_VERSION) return (FALSE);

                if (! read_long (fd, &pm1->stage, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->B_DONE, NULL)) goto readerr;
>>>>>>> Stashed changes
                if (! read_longlong (fd, &pm1->B, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->C_done, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->C_start, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->C, NULL)) goto readerr;
<<<<<<< Updated upstream
                /* "processed" is number of bits in stage 0, prime in stage 1 stored in pairs_done */
                if (! read_longlong (fd, &pm1->pairs_done, NULL)) goto readerr;
                if (! read_long (fd, &pm1->D, NULL)) goto readerr;
                if (! read_long (fd, &pm1->E, NULL)) goto readerr;
=======
                /* "processed" stored in bitarray_first_number */
                if (! read_longlong (fd, &pm1->bitarray_first_number, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->D, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->E, NULL)) goto readerr;
                if (! read_longlong (fd, &pm1->rels_done, NULL)) goto readerr;
>>>>>>> Stashed changes
                break;
        case 0x2c7330a8: // LL_MAGICNUM:
                //if (version != LL_VERSION) goto readerr;
                if (version != 1) goto readerr;

                // duplicated from commonb.c readLLSaveFile (minus reading data)
                w->work_type =  WORK_TEST;

                // Store data in pm1 as llhandle doesn't have two of the fields.
                // error_count is stored in E,
                // count (iterations) is stored in C,
                if (! read_long (fd, &pm1->E, NULL)) goto readerr;
                if (! read_long (fd, &pm1->C, NULL)) goto readerr;
                break;
        case 0x87f2a91b: // PRP_MAGICNUM:
                //if (version != PRP_VERSION) goto readerr;
                if (version != 4) goto readerr;

                w->work_type =  WORK_PRP;

                // Store data in pm1 as llhandle doesn't have two of the fields.
                // error_count is stored in E,
                // count (iterations) is stored in C,
                if (! read_long (fd, &pm1->E, NULL)) goto readerr;
                if (! read_long (fd, &pm1->C, NULL)) goto readerr;

                break;
        case 0x1567234D: // FACTOR_MAGICNUM:
                //if(version != FACTOR_MAGICNUM) goto readerr;
                if(version != 1) goto readerr;

                w->work_type =  WORK_FACTOR;
/* TODO implement. */

                break;
        default:
            goto readerr;
        }

/* Return success */
        _close (fd);
        return (TRUE);

readerr:
        _close (fd);
        return (FALSE);
}


int isTempFileName(char *filename)
{
/* Shortest reasonable file is p13_3 for 1*2^13+3 */
    if (strlen(filename) <= 4)
        return (FALSE);

/* Check if file is [mpef][0-9]*(_[0-9]*)(_[0-9]*) */
    char type = filename[0];
    if (!(type == 'm' || type == 'p' || type == 'e' || type == 'f'))
        return (FALSE);

    int underscores = 0;
    int i = 1;
    for (; i < strlen(filename); i++) {
        char d = filename[i];
        if (d == '_') {
            underscores++;
            if (underscores > 2)
                return (FALSE);
        } else if (d == '.' && i > 1) {
            break;
        } else if (d < '0' || d > '9') {
            return (FALSE);
        }
    }
    if (i == strlen(filename))
        return (TRUE);

    char * test = filename + i;

/* Check if file ends in optional .bu[0-9]* */
    if (strcmp(filename + i, ".bu") == 0) {
        i += 3;
        for (; i < strlen(filename); i++) {
            char d = filename[i];
            if (d < '0' || d > '9') {
                return (FALSE);
            }
        }
        return (TRUE);
    }

    return (FALSE);
}


int cmpFileName(void const *a, void const *b) {
    return strcmp((const char *) a, (const char *) b);
}


void restoreStatusMessage (
        char    *buf,
        unsigned int buflen)
{
        char    *orig_buf;

/* TODO augment by checking how many of these are mentioned in worktodo. */

        orig_buf = buf;

/* Get current working directory */
        char cwd[260];
        if (getcwd(cwd, sizeof(cwd)) == NULL) {
                sprintf(buf, BACKUP_CWD_ERROR);
                return;
        } else {
                char *last_dir = strrchr(cwd, '/');
                if (last_dir == NULL)
                    last_dir = cwd;
                else if (*last_dir == '/')
                    last_dir++;

                snprintf(buf, buflen, BACKUP_CWD, last_dir);
                buflen -= strlen(buf);
                buf += strlen(buf);
        }

/* Read up to MAX_BACKUP_FILES, 100 char backup filenames. */
        char status_filename[MAX_BACKUP_FILES][100] = {};

        struct dirent *dp;
        DIR *dfd;
        if ((dfd = opendir(".")) == NULL) {
                sprintf(buf, BACKUP_CWD_ERROR);
                return;
        }

        int status_files = 0;
        while ((dp = readdir(dfd)) != NULL)
        {
            if (dp->d_type != DT_REG)
                    continue;

            /* Assumes all backup filenames are < 100 characters . */
            if (strlen(dp->d_name) < 100 && isTempFileName(dp->d_name)) {
                    strcpy(status_filename[status_files++], dp->d_name);
                    if (status_files == MAX_BACKUP_FILES)
                            break;
            }
        }

/* Sort backup filenames. */
        qsort(status_filename, status_files, sizeof(status_filename[0]), cmpFileName);

        for (int i=0; i < status_files; i++) {
                char* filename = status_filename[i];

                // Load File into work_unit.
                struct work_unit w;
                pm1handle pm1;
                if (!restoreWorkUnitFromFile(filename, &w, &pm1)) {
                        snprintf(buf, buflen, BACKUP_PARSE_ERROR, filename);
                        buflen -= strlen(buf);
                        buf += strlen(buf);
                } else {

/* Process workunit and pm1 data into a status message */
                        char status[201] = {};
                        char *status_buf = status;

                        /* TODO should I print the number here? And why is K a float? */
                        // status_buf += sprintf(status_buf, "%g*%lu^%lu%c%lu | ", w.k, w.b, w.n, w.c < 0 ? '-' : '+', abs(w.c));
                        switch (w.work_type) {
                        case WORK_ECM:
                                // TODO print bounds?
                                status_buf += sprintf(status_buf, "ECM | Curve %d | Stage %ld (%.1f%%)",
                                    w.curves_to_do, pm1.stage + 1, 100 * w.pct_complete);
                                break;
                        case WORK_PMINUS1:
                                switch (pm1.stage) {
                                case 3: //PM1_STAGE3
/* Stage 1, pairs_done = processed = bit_number */
                                    status_buf += sprintf(status_buf, "P-1 | Stage 1 (%.1f%%) B1 <%lu",
                                        100 * w.pct_complete, pm1.pairs_done);
                                    break;

                                case 0: //PM1_STAGE0
/* Stage 1 after small primes, pairs_done = processed = prime */
                                    status_buf += sprintf(status_buf, "P-1 | Stage 1 (%.1f%%) B1 @ %lu",
                                        100 * w.pct_complete, pm1.pairs_done);
                                    break;

                                case 1: //PM1_STAGE1
/* Stage 2 after small primes, pairs_done = processed = B1 bound */
                                    status_buf += sprintf(status_buf, "P-1 | B1=%0.f complete, Stage 2 (%.1f%%)",
                                        (double) pm1.B, 100 * w.pct_complete);
                                    break;

                                case 2: //PM1_DONE
/* P-1 done */
                                    status_buf += sprintf(status_buf, "P-1 | B1=%.0f", (double) pm1.B);
                                    if (pm1.C > pm1.B) {
                                        status_buf += sprintf(status_buf, ",B2=%.0f", (double) pm1.C);
                                        if (pm1.E >= 2)
                                            status_buf += sprintf(status_buf, ",E=%lu", pm1.E);
                                    }
                                    status_buf += sprintf(status_buf, " complete");
                                    break;
                                }
                                break;

                        case WORK_TEST:
                                status_buf += sprintf(status_buf, "LL  | Iteration %lu/%lu [%0.2f%%]",
                                    pm1.C, w.n, 100 * w.pct_complete);
                                break;

                        case WORK_PRP:
                                status_buf += sprintf(status_buf, "PRP | Iteration %lu/%lu [%0.2f%%]",
                                    pm1.C, w.n, 100 * w.pct_complete);
                                break;

                        case WORK_FACTOR:
                                break;

                        }

                        if (strlen(status) == 0)
                            sprintf(status, "UNKNOWN");

                        snprintf(buf, buflen, BACKUP_STATUS, filename, status);
                        buflen -= strlen(buf);
                        buf += strlen(buf);
                }
        }
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
