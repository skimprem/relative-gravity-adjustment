from os import path
import sys
import argparse
import tkinter as tk
from tkinter import filedialog as fd
from tkinter import simpledialog as sd
from tkinter import messagebox as mb
from tkinter import scrolledtext as st
from src.dataloader import DataLoader
import pandas as pd

# GUI = False
GUI = True 

if sys.platform.startswith('win32'):
    GUI = True

if GUI:
    data_file = fd.askopenfilename(
        defaultextension='.dat',
        filetypes=[('CG-x data files', '*.dat'), ('All files', '*')],
        title='Choose data file'
    )
    degree = sd.askinteger(
        title='Degree input',
        prompt='Input degree',
        initialvalue=1,
        minvalue=1,
        maxvalue=6
    )
    anchor = sd.askstring(
        title='Anchor input',
        prompt='Input anchor'
    )
    plot = mb.askyesno(
        title='Plotting',
        message='Want to make a plot?'
    )
    verbose = mb.askyesno(
        title='Report',
        message='Want to print a report?'
    )
else:
    parser=argparse.ArgumentParser(
        prog='rga',
        description='Read CG-6 data file and compute ties',
        epilog='This program read CG-6 data file, then compute ties by ...',
        exit_on_error=False
    )

    parser.add_argument('data_file', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('--meter', metavar='meter_type', help='Type of meter')
    parser.add_argument('--degree', metavar='degree', help='Degree')
    parser.add_argument('--verbose', action='store_true', help='Verbose mode')
    parser.add_argument('--anchor', metavar='anchor', help='Anchor')
    parser.add_argument('--plot', action='store_true', help='Get plot')

    args=parser.parse_args()
    data_file = args.data_file[0].name
    degree = int(args.degree)
    if args.anchor:
        anchor = args.anchor
    else:
        anchor = None
    if args.verbose:
        verbose = args.verbose
    else:
        verbose = None
    if args.plot:
        plot = args.plot
    else:
        plot = None

data = DataLoader.load('CG6', data_file)
if anchor:
    result = data.invert(degree=degree, anchor=anchor)
else:
    result = data.invert(degree=degree)

anchor = result.differences['anchor']

report_table = pd.DataFrame(columns=['From', 'To', 'Tie (uGals)', 'Err (uGals)'])

for difference in result.differences['differences']:
    report_table = pd.concat([report_table, pd.DataFrame([[anchor, difference[0], difference[1], difference[2]]], columns=['From', 'To', 'Tie (uGals)', 'Err (uGals)'])])

report = report_table.to_markdown(
    index=False,
    headers=['From', 'To', 'Tie (uGals)', 'Err (uGals)'],
    tablefmt='simple',
    floatfmt='.2f'
)

basename = path.splitext(data_file)[0].split('/')[-1]

if verbose:
    if GUI:
        # Creating tkinter window
        win = tk.Tk()
        win.title('Report')
          
        # Creating scrolled text area
        # widget with Read only by
        # disabling the state
        text_area = st.ScrolledText(win,
                                    width = 50, 
                                    height = 10
        )
          
        text_area.grid(column = 0, pady = 10, padx = 10)
          
        # Inserting Text which is read only
        text_area.insert(tk.END, report)
          
        # Making the text read only
        text_area.configure(state ='disabled')
        win.mainloop()
    else:
        print(report)

default_output_file_report = 'report_'+basename+'.txt'

if GUI:
    output_file_report = fd.asksaveasfilename(
        defaultextension=".txt",
        filetypes=[('ACSII text file', '*.txt'), ('All files', '*')],
        initialfile=default_output_file_report,
        title="Save Report")
else:
    output_file_report = default_output_file_report

with open(output_file_report, 'w', encoding='utf-8') as report_file:
    report_file.write(report)
    report_file.close()

if plot:
    if GUI:
        output_plot = fd.asksaveasfilename(
            defaultextension='.pdf',
            filetypes=[('pdf', '*.pdf'), ('All files', '*')],
            initialfile=basename+'.pdf',
            title='Save Plot'
        )
    else:
        output_plot = basename+'.pdf'
result.plot(output_plot)