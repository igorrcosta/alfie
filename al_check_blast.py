from sys import argv

''' python al_check_blast *blast_output.m8* *total number of ALs*
Prints the ids of the ALs that were not found in
a Blast search with -m8 output setting''' 

def check_blast(blast_m8, size):
    with open(blast_m8, 'r') as blast:
        list_n = range(1, size)
        for l in blast:
            n = int(l.split('|')[0].split('_')[1])
            if n in list_n:
                list_n.remove(n)
    print list_n

if __name__ == '__main__':
    blast_m8 = argv[1]
    if len(argv) > 2 and type(argv[2]) == int:
        size = argv[2]
    else:
        size = 1000
    check_blast(blast_m8, size)
