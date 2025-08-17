import os
from tabulate import tabulate


def print_with_decoration(info, decoration_char="*", t_width_max=80):
    try:
        t_width, t_height = os.get_terminal_size()
    except OSError:
        t_width, t_height = 80, 24
    t_width = min(t_width, t_width_max)

    if "\n" in info:
        info_list = info.split("\n")
    else:
        info_list = [str(info)]
    for info in info_list:
        if len(info) > t_width - 10:
            info_to_print = info[: t_width - 10] + "..."
        else:
            info_to_print = info
        decoration_len = (t_width - len(info_to_print) - 2) // 2
        print(
            f"{decoration_char * decoration_len} {info_to_print} {decoration_char * decoration_len}"
        )


def print_decoration_line(decoration_char="*"):
    try:
        t_width, t_height = os.get_terminal_size()
    except OSError:
        t_width, t_height = 80, 24
    t_width = min(t_width, t_width_max)
    print(decoration_char * t_width)


def str_mappable(mappable, **kwargs):
    table = []
    for key, value in mappable.items():
        table.append([key, str(value)])
    return str(tabulate(table, **kwargs))


def print_mappable(mappable, **kwargs):
    table = []
    for key, value in mappable.items():
        table.append([key, str(value)])
    print(tabulate(table, **kwargs))
