import csv
from collections import defaultdict
dyes = {}
name_count = defaultdict(int)
wavelength =[]
fname = 'McNamara_Spectra.csv'

with open(fname, 'r') as infile:
    spectra_reader = csv.reader(infile)
    row_count = 0
    row = spectra_reader.next()
    while row[0]!='dye name':
        row=spectra_reader.next()
        print row[0]

    dye_names  = row[1:]
    nice_dye_names = []
    columns_to_dyename = {}
    for ci, dye in enumerate(dye_names):
        nice_name = "_".join(dye.split())
        nice_name = "-".join(nice_name.split('/'))
        nice_dye_names.append(nice_name)
        name_count[nice_name]+=1
        if name_count[nice_name]==1:
            columns_to_dyename[ci] = (nice_name, 0)
        elif name_count[nice_name]==2:
            columns_to_dyename[ci] = (nice_name, 1)
        else:
            columns_to_dyename[ci] = ('', -1)

            
        dyes[nice_name]=[[], []]

    while row[0]!='57':
        row=spectra_reader.next()


    for row in spectra_reader:
        wavelength.append(int(row[0]))
        for ci, val in enumerate(row[1:]):
            dye_name, exc_emn = columns_to_dyename[ci] 
            if exc_emn>0:
                try:
                    dyes[dye_name][exc_emn].append(float(val))
                except:
                    dyes[dye_name][exc_emn].append(0.0)


for dye_name in dyes:
    if len(dyes[dye_name][0]):
        with open('SpectraLibrary/'+dye_name+'.exc', 'w') as outfile:
            if len(wavelength)==len(dyes[dye_name][0]):
                for wl, exc in zip(wavelength, dyes[dye_name][0]):
                    if exc>0:
                        outfile.write(str(wl)+'\t'+str(exc)+'\n')
            else:
                print dye_name, "spectrum does not have the right length"

    if len(dyes[dye_name][1]):
        with open('SpectraLibrary/'+dye_name+'.emn', 'w') as outfile:
            if len(wavelength)==len(dyes[dye_name][1]):
                for wl, emn in zip(wavelength, dyes[dye_name][1]):
                    if emn>0:
                        outfile.write(str(wl)+'\t'+str(emn)+'\n')
            else:
                print dye_name, "spectrum does not have the right length"

