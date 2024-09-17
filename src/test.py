import csv

with open("src/tfbs.csv", "r") as csvfile:
    datareader = csv.reader(csvfile)
    header = next(datareader)  # yield the header row
    count = 0
    for row in datareader:
        print(row)