from os import path
import argparse
from src.dataloader import DataLoader
import pandas as pd

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
parser.add_argument('--meter_type', metavar='meter_type', help='Meter type')
parser.add_argument('--plot', action='store_true', help='Get plot')

args=parser.parse_args()

data_file = args.data_file[0].name
degree = int(args.degree)

data = DataLoader.load(args.meter_type, data_file)
if args.anchor:
    result = data.invert(degree=degree, anchor=args.anchor)
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

with open('report_'+basename+'.txt', 'w', encoding='utf-8') as report_file:
    report_file.write(report)
    report_file.close()

if args.verbose:
    print(report)

if args.plot:
    result.plot(basename+'.pdf')