import numpy as np
import os
import gzip
import re
import subprocess


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def gzipped_extension(filename):
    """
    Return true if file seems to be gzipped based on extension
    """
    return filename[-3:] == '.gz'

def get_uncompressed_filename(filename):
    """
    For .gz file, returns the filename when the file is uncompressed.
    """
    return re.sub('[.]gz$', '', filename)

def gzip_file(filename):
    if which('pigz') is not None:
        subprocess.call('pigz -f %s' % filename, shell=True)
    else:
        subprocess.call('gzip %s' % filename, shell=True)

def open_file(file_name, mode=None):
    if gzipped_extension(file_name):
        if mode is None:
            mode = 'rt'
        return gzip.open(file_name, mode)
    else:
        if mode is None:
            mode = 'r'
        return open(file_name, mode)

def write_list(l, filename):
    with open(filename, 'w') as output_file:
        for i in l:
            output_file.write('%s\n' % i)

def read_list(filename):
    return np.array([x.strip() for x in open(filename)])

def _get_aux_files(mtx_file):
	base_name = re.sub('[.]gz$', '', mtx_file)
	base_name = re.sub('[.]mtx$', '', base_name)

	rows_file = '%s.rows.txt' % base_name
	columns_file = '%s.columns.txt' % base_name

	return (rows_file, columns_file)

def write_mtx_file(mat, row_names, column_names, filename):
    from scipy.io import mmwrite

    rows_file, columns_file = _get_aux_files(filename)
    write_list(list(row_names), rows_file)
    write_list(list(column_names), columns_file)
    
    if filename[-3:] == '.gz':
        uncompressed_file = get_uncompressed_filename(filename)
        mmwrite(target=uncompressed_file, a=mat)
        gzip_file(uncompressed_file)
    else:
        mmwrite(target=filename, a=mat)

def load_mtx_file(mtx_file):
    from scipy.io import mmread

    mat = mmread(mtx_file)
    features_file, cells_file = _get_aux_files(mtx_file)

    if os.path.exists(features_file):
        raise ValueError('%s features file not found when loading %s mtx file.' % (features_file, mtx_file))
    if os.path.exists(cells_file):
        raise ValueError('%s cells file not found when loading %s mtx file.' % (cells_file, mtx_file))

    return (mat, read_list(features_file), read_list(cells_file))

