"""Inspired by http://stackoverflow.com/a/34325723/6600524"""
import sys

lastPercents = -1.

# Print iterations progress
def printProgress (iteration, total, prefix = 'Progress:', suffix = 'Complete', decimals = 1, barLength = 50):
    global lastPercents
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    formatStr       = "{0:." + str(decimals) + "f}"
    percents        = formatStr.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    if percents != lastPercents:
        sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
        lastPercents = percents
    sys.stdout.flush()
    if iteration == total:
        lastPercents = -1.
        sys.stdout.write('\n')
        sys.stdout.flush()