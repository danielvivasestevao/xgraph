import os
import re


def file_exists(path_to_file: str):
    """
    :param path_to_file: path to a file
    :return: True if the file exists, False otherwise
    """
    return os.path.isfile(path_to_file)


def ensure_dir(path):
    """
    Creates a directory if it does not exist already.

    :param path: path to a directory
    """
    # ensure directory exists
    d = os.path.dirname(path + '/')
    if not os.path.exists(d):
        os.makedirs(d)


def get_safe_filename(filename: str):
    """
    Checks if a file exists, and if it does, increments the last number in its
    name and returns the full path with the new name. Returns the original path
    and filename otherwise.

    :param filename: path to a potential file
    :return: the given path if it does not exist; an increment of its filename
     including the path otherwise
    """
    if os.path.isfile(filename):
        return get_safe_filename(increment(filename))
    else:
        return filename


# https://code.activestate.com/recipes/442460-increment-numbers-in-a-string/
def increment(file_name: str):
    last_num = re.compile(r'(?:[^\d]*(\d+)[^\d]*)+')
    m = last_num.search(file_name)
    if m:
        nxt = str(int(m.group(1))+1)
        start, end = m.span(1)
        file_name =\
            file_name[:max(end-len(nxt), start)] + nxt + file_name[end:]
    return file_name
