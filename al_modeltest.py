from subprocess import Popen, PIPE
import os
import shlex
import argparse

def argument_parser():

    parser = argparse.ArgumentParser(description = 'al_modeltest.py',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-o', '--oupath', nargs = '?', type = str, default = 'models.txt',\
                        dest = 'outpath', help = 'Path where the modeltest results will be saved. (default: %(default)s)')
    parser.add_argument('-i', '--inpath', nargs = '?', type = str, required = True,\
                        dest = 'inpath', help = 'Path with the sample.log and nexus files.')
    args = parser.parse_args().__dict__
    return args

def main(args):
    if args['inpath'][-1] != '/':
	args['inpath'] += '/'
    model = {}
    table = ''
    with open(args['inpath'] + 'sample.log', 'r') as sample:
        file_list = sample.readlines()
    total = len(file_list)
    n = 0
    for f in sorted(file_list):
        f = f[:-1]
        print('Runing modeltest: ', f, str(n+1), '/', str(total))
        n += 1
        command = 'java -jar jModelTest.jar -d ' + args['inpath'] + f + ' -f -i -g 4 -s 11 -BIC -a' #-DT -AIC -AICc
        a = com(command)
        result = a[0][-1].split() #test model a c g t kappa titv
        try:
            model[result[1]] += 1
	except:
            model[result[1]] = 1
        table += '\t'.join([f.split('.')[0], result[7], result[1]]) + '\n'
        with open(args['inpath'] + 'modeltest.'+f, 'w') as outfile:
            for l in a[0]: 
                outfile.write(l + '\n')
    with open(args['outpath'], 'w') as models:
        for k in model.keys():
            models.write(k +' '+ str(model[k]) + '\n')
        models.write(table)

def com(command):
    o = Popen(command, shell = True, stdout=PIPE, stderr = PIPE)
    with o.stdout as out:
        r = [l[:-1] for l in out]
    with o.stderr as err:
        e = [l[:-1] for l in err]
    return r, e

if __name__ == '__main__':
    main(argument_parser())
