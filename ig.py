#!/usr/bin/env python

#input generator
#takes parameters from file argument and creates a bunch of input files!
import subprocess
import os
import re
import time
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', nargs=1)
parser.add_argument('--n_max_sims', nargs=1, default='3')
parser.add_argument('--use_log_kappa', nargs=1, default='False')
args = parser.parse_args(sys.argv[1:])

input_name = args.input_file[0]

n_max_sims = int(args.n_max_sims[0])
use_log_kappa = False
if (args.use_log_kappa[0] == 'True'):
  use_log_kappa = True

#input_name = sys.argv[1]
print(input_name)


required_parameters = ['kappas', 'etas', 'Ls']
auto_parameters = ['max_states', 'total_states']

def make_dir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
#input etc
def general_input_converter(input_name):
    result = {}
    f = open(input_name)
    lines = f.readlines()
    f.close()
    for l in lines:
        l = l.strip()
        l = re.sub('=', ' ', l)
        l = re.sub(' +', ' ', l)
        z = l.split(' ')
        result[z[0]] = z[1:]
    return result

p = general_input_converter(input_name)
    
def convert_to_number(a):
    offset = 0
    if a[0] == 'o': 
        offset = 1
        return (1.0 / float(a[offset:]))
    b = a[offset:]
    if '.' in b:
        if int(b.split('.')[-1]) == 0:
            return int(b.split('.')[0])
        else:
            return float(b)
    else:
        return int(b)

def determine_parameters(kappa, eta, L):
    result = {}
    kn = convert_to_number(kappa)
    hn = convert_to_number(eta)

    max_states = "80"
    #max_states = "5"
    if (kn != 1 or hn != 1):
        max_states += " 150"
        #max_states += " 3"

    temp = 250
    if (hn != 1):
        if (eta > 80):
            temp += 75
        elif (eta > 200):
            temp += 150
        else:
            temp += 0
        #max_states += " 3"
    if (kn > 12.0):
        temp += 50
    if (kn > 13.0):
        temp += 50

    max_states += ' ' + str(temp)

    result['max_states'] = max_states

    total_states = 0
    for m in max_states.split():
        total_states += int(m)
    result['total_states'] = str(total_states)

    try:
        result["sf_channels"] = p["sf_channels"]
    except KeyError:
        1

    return result


name_prefix = "test1"
sims = {}
for kappa in p['kappas']:
    for eta in p['etas']:
        for L in p['Ls']:
            a = determine_parameters(kappa, eta, L)
            name = kappa + '_' + eta + '_' + L;

            content = ""
            content += "mass_ratios=" + str(convert_to_number(kappa)) + ';\n'
            content += "omega_ratios=" + str(convert_to_number(eta)) + ';\n'
            content += "ltotal=" + str(convert_to_number(L)) + ';\n'
            for auto_p in auto_parameters:
                content += auto_p + "=" + a[auto_p] + ';\n'


            sims[name] = {}
            for par in a:
                if (par not in auto_parameters):
                    if (par not in required_parameters):
                        content += par + "=" + " ".join(a[par]) + ';\n'
                        sims[name][par] = a[par]


            sims[name]['content'] = content

for skey in sims:
    s = sims[skey]
    d = input_name.split('.')[0] + "_dir/"
    make_dir(d)
    s['dir'] = d
    s['file_main'] = d + skey
    w = open(s['file_main'], 'w+')
    w.write(s['content'])



for skey in sims:
    s = sims[skey]
    s['file_states'] = s['dir'] + skey + ".states"
    s['file_out'] = s['dir'] + skey + ".out"
    s['file_sf'] = s['dir'] + skey + ".sf"
    s['file_info'] = s['dir'] + skey + ".info"
    s['finished'] = False
    s['active'] = False
    if os.path.isfile(s['file_states']):
        f = open(s['file_states'])
        lines = f.readlines()
        f.close()
        not_finished = False
        if not os.path.isfile(s['file_info']):
            not_finished = True

        for chan in s['sf_channels']:
            if not os.path.isfile(s['file_sf'] + chan):
                not_finished = True

        s['finished'] = not not_finished
            


n_active_sims = 0
i = 0
active_sims = []
while (i < n_max_sims):
    active_sims.append({})
    i += 1


def fill_active_sims(active_sims):
    flag = False
    i = 0
    while (i < n_max_sims):
        a = active_sims[i]
        if (a == {}):
            flag = activate_next_sim(active_sims, i)
        elif (a['proc'].poll() != None):
            a['out_stream'].close()
            sims[a['key']]['finished'] = True
            sims[a['key']]['active'] = False
            flag = activate_next_sim(active_sims, i)
            f = open(sims[a['key']]['file_info'], 'w+')
            f.write(str(sims[a['key']]))

        i += 1
        if (flag == True): break

    return flag

def activate_next_sim(active_sims, i):
    a = active_sims[i]
    all_finished = True
    for skey in sims:
        s = sims[skey]

        if s['finished'] == False:
            all_finished = False
            if s['active'] == False:
                s['active'] = True
                a['out_stream'] = open(s['file_out'], 'w+')
                command = ['./3imb', s['file_main']]
                a['proc'] = subprocess.Popen(command, stdout=a['out_stream'])
                a['key'] = skey
                return False

    return all_finished

            

def report_active_sims(active_sims):
    for a in active_sims:
        if (a == {}):
            print "EMPTY"
        else:
            output = a['key'] + " " + str(sims[a['key']]['active'])
            print output
    
    print " "

timer = 0
while (1):
    if (not fill_active_sims(active_sims)):
        report_active_sims(active_sims)
        time.sleep(10)
        timer += 10
        print "t = " + str(timer)
    else:
        print "Finished! Timer = " + str(timer)
        break

graphs_sf_out = open(s['dir'] + "graphs_sf.sh", 'w+')
for chan in p['sf_channels']:
    command = ['python', 'g3.py']
    for skey in sims:
        s = sims[skey]
        command.append(s['file_main'])

    command.append('--suffix')
    command.append('.sf' + chan)
    command.append('--output')
    command.append(s['dir'] + 'graph' + '_sf' + chan + '.pdf')

    f = open('/dev/null')
    subprocess.call(command, stdout=f)
    f.close()
    print chan


