import h5py
import numpy as np

h5 = h5py.File('even1k.weighted.dm.pc')

# load results from "permanova_sample_filtering.R"
print('Reading keepers')
with open('permanova_to_keep.txt', 'r') as infile:
    infile.readline()  # skip header
    keepers = [x[:-1].encode('utf-8') for x in infile]

indices = []
skipped = []
print('Finding indices')
biglist = h5['order'][0:]
for i,x in enumerate(keepers):
    check = np.where(biglist==x)
    if len(check) > 0 and len(check[0]) > 0:
        indices.append(check[0][0])
    else:
        skipped.append(x)
print(f'Found {len(indices)} indices, missed {len(skipped)}')
print('writing')
with open('keeper_indices.txt','w') as outfile:
    towrite = [str(x) for x in indices]
    outfile.write('\n'.join(towrite))
    outfile.write('\n')

# with open('keeper_indices.txt','r') as outfile:
#     indices = [int(x[:-1]) for x in outfile]

indices.sort()

with h5py.File('even1k.weighted.dm.pc','r') as h5:
    biglist = h5['order'][0:]

filtered_names = biglist[indices]
with open('keeper_names_sorted.txt','w') as outfile:
    towrite = [x.decode('utf-8') for x in filtered_names]
    outfile.write('\n'.join(towrite))
    outfile.write('\n')

with h5py.File('even1k.weighted.dm.pc','r') as h5:
    print('filtering A')
    with h5py.File('half_filtered.weighted.dm.pc','w') as outfile:
        outfile.create_dataset('matrix', dtype=np.single, data=h5['matrix'][:, indices])
    print("done the first half!")

with h5py.File('half_filtered.weighted.dm.pc','r') as h5:
    print('filtering B')
    with h5py.File('filtered.weighted.dm.pc','w') as outfile:
        outfile.create_dataset('matrix', dtype=np.single, data=h5['matrix'][indices, :])

    print("DONE part 2!!")
