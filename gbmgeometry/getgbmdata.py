from __future__ import print_function

import datetime
import ftplib
import sys
import time


class GetGBMData(object):
    def __init__(self, trigger=None, daily=None):

        self.ftp = ftplib.FTP('legacy.gsfc.nasa.gov', 'anonymous', 'crap@gmail.com')

        if trigger is not None:

            self._type = 'triggered'

            year = '20' + trigger[:2]
            self._directory = 'fermi/data/gbm/triggers/' + year + '/bn' + trigger + '/current'

        elif daily is not None:

            self._type = 'daily'

            date = daily.split('/')

            self._directory = 'fermi/data/gbm/daily/20' + date[2] + '/' + date[0] + '/' + date[1] + '/current'


        else:

            print("You did something wrong!")
            return

        try:
            self.ftp.cwd(self._directory)
        except ftplib.error_perm:
            print(self._directory)
            print("Awww snap! This data entry does not exist at the FSSC. Exiting!\n")
            self.ftp.quit()
            return

        self._file_list = self.ftp.nlst()

        self._detectors = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']

        self._where = ''

    def select_detectors(self, *dets):
        self._detectors = dets

        self._num_dets = len(dets)

    def set_destination(self, destination):

        self._where = destination

    def _get(self, items):

        n_items = len(items)

        progress = ProgressBar(n_items)

        progress_bar_iter = max(int(n_items / 100), 1)

        for i, item in enumerate(items):

            if i % progress_bar_iter == 0:
                progress.animate((i + 1))

            self.ftp.retrbinary('RETR ' + item, open(self._where + item, 'wb').write)

        progress.animate(n_items)

    def get_rsp_cspec(self):

        to_get = []

        for i in self._file_list:
            for j in self._detectors:
                if ".rsp" in i and "cspec" in i and j + '_' in i:
                    to_get.append(i)
        self._get(to_get)

    def get_rsp_ctime(self):

        to_get = []

        for i in self._file_list:
            for j in self._detectors:
                if ".rsp" in i and "ctime" in i and j + '_' in i:
                    to_get.append(i)
        self._get(to_get)

    def get_ctime(self):

        to_get = []

        for i in self._file_list:
            for j in self._detectors:
                if "ctime" in i and j + '_' in i and 'rsp' not in i:
                    to_get.append(i)
        self._get(to_get)

    def get_cspec(self):

        to_get = []

        for i in self._file_list:
            for j in self._detectors:
                if "cspec" in i and j + '_' in i and 'rsp' not in i:
                    to_get.append(i)

        self._get(to_get)

    def get_trigdat(self):

        to_get = []

        for i in self._file_list:
            if "trigdat" in i:
                to_get.append(i)

        self._get(to_get)

    def get_tte(self):

        to_get = []

        for i in self._file_list:
            for j in self._detectors:
                if "tte" in i and j + '_' in i:
                    to_get.append(i)

        self._get(to_get)

    def get_tcat(self):

        to_get = []

        for i in self._file_list:
            if "tcat" in i:
                to_get.append(i)

        self._get(to_get)

    def get_plots(self):

        to_get = []

        for i in self._file_list:
            if ".gif" in i or ".pdf" in i:
                to_get.append(i)

        self._get(to_get)

    def get_poshist(self):

        to_get = []

        for i in self._file_list:
            if "poshist" in i:
                to_get.append(i)

        self._get(to_get)

    def get_spechist(self, det='g'):

        to_get = []

        for i in self._file_list:
            for j in det:
                if "spechist" in i and '_' + j + '_' in i:
                    to_get.append(i)

        self._get(to_get)


class ProgressBar:
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 50
        self.startTime = time.time()
        self.lastIter = 0
        self.__update_amount(0)

    def animate(self, iter):

        try:

            print('\r', self, end='')
            sys.stdout.flush()
            self.lastIter = iter
            self.update_iteration(iter + 1)

        except:
            # Do not crash in any case. This isn't an important operation
            pass

    def increase(self):

        self.animate(self.lastIter + 1)

    def _check_remaining_time(self, delta_t):

        # Seconds per iterations
        s_per_iter = delta_t / float(self.lastIter)

        # Seconds to go (estimate)
        s_to_go = s_per_iter * (self.iterations - self.lastIter)

        # I cast to int so it won't show decimal seconds

        return str(datetime.timedelta(seconds=int(s_to_go)))

    def update_iteration(self, elapsed_iter):

        delta_t = time.time() - self.startTime

        elapsed_iter = min(elapsed_iter, self.iterations)

        if elapsed_iter < self.iterations:

            self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
            self.prog_bar += '  %d / %s in %.1f s' % (elapsed_iter, self.iterations, delta_t)
            self.prog_bar += ' (%s remaining)' % self._check_remaining_time(delta_t)

        else:

            self.__update_amount(100)
            self.prog_bar += '  completed in %.1f s' % (time.time() - self.startTime)

    def __update_amount(self, new_amount):
        percent_done = min(int(round((new_amount / 100.0) * 100.0)), 100)
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
                        (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)
