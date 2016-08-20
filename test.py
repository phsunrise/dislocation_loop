import sys

with open("test.txt", 'w') as f:
    f.write(str(len(sys.argv)))
    f.write('\n')
    f.write('  '.join(sys.argv))
