import curses
from curses import wrapper
from curses.textpad import Textbox, rectangle
import pandas as pd
import time
import pprint

#stdscr = curses.initscr

def generate_set(n_atoms, idof, bonds, angles, linear_angles, out_of_plane, dihedrals) -> dict:
    return 0

def write_generated_IC(bonds, angles, linear_angles, out_of_plane, dihedrals, idof) -> str:
    return ""

def main(stdscr):
    curses.init_pair(1, curses.COLOR_BLUE, curses.COLOR_YELLOW)
    curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED, curses.COLOR_WHITE)
    BLUE_AND_YELLOW = curses.color_pair(1)
    GREEN_AND_BLACK = curses.color_pair(2)
    ORANGE_AND_WHITE = curses.color_pair(3)
    stdscr.nodelay(True)


    x,y = 0,0
    while True:
        try:
            key = stdscr.getkey()
        except:
            key = None
        if key == "KEY_LEFT":
            x -= 1
        elif key == "KEY_RIGHT":
            x += 1
        elif key == "KEY_UP":
            y -= 1
        elif key == "KEY_DOWN":
            y += 1
        
        stdscr.clear()
        stdscr.addstr(y,x, "Segmentation fault, you loser")
        stdscr.refresh()

wrapper(main)