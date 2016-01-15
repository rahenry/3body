#!/usr/bin/env python

#input generator
#takes parameters from file argument and creates a bunch of input files!
import subprocess
import os
import re
import time
import sys
import math

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', nargs=1)
parser.add_argument('--make_graphs', nargs='?', default=None, const='1')
parser.add_argument('--run_project', nargs='?', default=None, const='1')
parser.add_argument('--max_active', nargs='?', default='4', const='4')
parser.add_argument('--finish_all', nargs='?', default=None, const='1')

args = parser.parse_args(sys.argv[1:])

input_name = os.path.abspath(args.input_file[0])

#input_name = sys.argv[1]


required_parameters = ['kappas', 'etas', 'Ls']
auto_parameters = ['max_states', 'total_states']

def make_dir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    return path

directory = input_name.split('.')[0] + "_dir/"
make_dir(directory)
dir_3body = ""
for x in directory.split('/')[1:-3]:
    dir_3body += '/' + x
dir_3body += '/'

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

def determine_parameters(s):
    result = {}
    kn = convert_to_number(s['kappa'])
    hn = convert_to_number(s['eta'])

    current_states = 0
    if os.path.isfile(s['file_states']):
      f = open(s['file_states'], 'rw+')
      current_states = len(f.readlines())
      f.close()
    max_states = "80"
    total_states = 0
    total_states += int((max_states.split())[-1])
    #max_states = "5"
    if (kn != 1 or hn != 1):
        max_states += " 150"
        total_states += int((max_states.split())[-1])
        #max_states += " 3"
    if (hn != 1):
        if (args.finish_all and current_states > 50):
          max_states += ' ' + str(current_states - total_states)
        else:
          if (s['eta'] > 240):
            max_states += " 550"
          elif (s['eta'] > 160):
            max_states += " 450"
          elif (s['eta'] > 80):
            max_states += " 350"
          else:
            max_states += " 250"

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
            name = kappa + '_' + eta + '_' + L;
            sims[name] = {}
            sims[name]['eta'] = eta
            sims[name]['kappa'] = kappa
            sims[name]['L'] = L


for skey in sims:
    s = sims[skey]
    d = directory
    s['dir'] = d
    s['file_main'] = d + skey
    s['file_states'] = d + skey + ".states"
    s['file_out'] = d + skey + ".out"
    s['file_sf'] = d + skey + ".sf"
    s['file_info'] = d + skey + ".info"
    s['finished'] = True
    s['active'] = False
    s['file_job'] = d + skey + ".job"


    a = determine_parameters(s)
    kappa = s['kappa']
    eta = s['eta']

    content = ""
    content += "mass_ratios=" + str(convert_to_number(kappa)) + ';\n'
    content += "omega_ratios=" + str(convert_to_number(eta)) + ';\n'
    content += "ltotal=" + str(convert_to_number(L)) + ';\n'
    for auto_p in auto_parameters:
        content += auto_p + "=" + a[auto_p] + ';\n'

    for par in a:
        if (par not in auto_parameters):
            if (par not in required_parameters):
                content += par + "=" + " ".join(a[par]) + ';\n'
                s[par] = a[par]

    s['content'] = content
    w = open(s['file_main'], 'w+')
    w.write(s['content'])
    #if (not s['finished']): print s['file_out']

    if os.path.isfile(s['file_states']):
        f2 = open(s['file_states'])
        f2l = f2.readlines()
        if (f2l < a['max_states']):
            s['finished'] = False
    else:
        s['finished'] = False


    for chan in s['sf_channels']:
        if not os.path.isfile(s['file_sf'] + chan):
            s['finished'] = False
        

f = open("job_base.pbs")
content_job_base = ""
for l in f.readlines():
    content_job_base += l


for skey in sims:
    s = sims[skey]
    s['content_job'] = content_job_base
    s['content_job'] += "cd " + s['dir'] + '\n'
    s['content_job'] += "../../3imb " + s['file_main'] + ' > ' + s['file_out'] + '\n'
    f = open(s['file_job'], 'w+')
    f.write(s['content_job'])

start_script = open(directory + "start_qsub.sh", 'w+')
for skey in sims:
    s = sims[skey]
    if not s['finished']:
       command = 'qsub ' + s['file_job']
       start_script.write(command + "\n")
start_script.close()

f = open(input_name)
kappa_log = False
eta_log = False
lines = f.readlines()
if "kappa_log" in lines[0]:
    kappa_log = True
if "eta_log" in lines[0]:
    eta_log = True
f.close()
  

def make_graphs():
    print("Generating graphs...")
    graphs_sf = open(directory + "graphs_sf.sh", 'w+')
    for chan in p['sf_channels']:
        command = "python " + dir_3body + "g3.py "
        for skey in sims:
            s = sims[skey]
            if os.path.isfile(s['file_sf'] + chan):
                command += s['file_main'] + ' '

        command += '--suffix ' + '.sf' + chan + ' ' + '--output ' + s['dir'] + 'graph_sf' + chan + '.pdf'
        graphs_sf.write(command + '\n')
    graphs_sf.close()
        
    # energy graph
    command = "python " + dir_3body + "g3.py "
    energy_files = generate_energy_data()

def get_key(thing):
    return thing[0]
    
def generate_energy_data():
    energy_dir = make_dir(d + "energy_data/")
    energy_graph = {}
    files = {}

    for eta in p['etas']:
        energy_graph[eta] = []
        for skey in sims:
            s = sims[skey]
            if (not s['finished']): continue
            if (eta == s['eta']):
                file_name = s['file_main'] + ".energies"
                f = open(file_name)
                data = (f.readlines()[-1].split())
                energy = float(data[1])
                f.close()
                kappa = float(s['kappa'])
                if (kappa_log): kappa = math.log(kappa)
                energy_graph[eta].append((kappa, energy))
        energy_graph[eta].sort(key=get_key)

    files['etas'] = []
                
    for etakey in energy_graph:
        etagraph = energy_graph[etakey]
        #print(etagraph)
        file_name = energy_dir + "eta_" + etakey + ".energies"
        f = open(file_name, 'w+')
        files['etas'].append(file_name)
        for graph_point in etagraph:
            kappa = graph_point[0]
            energy = graph_point[1]
            f.write(str(kappa))
            f.write(" ")
            f.write(str(energy))
            f.write('\n')
        f.close()

    energy_graph = {}

    for kappa in p['kappas']:
        energy_graph[kappa] = []
        for skey in sims:
            s = sims[skey]
            if (not s['finished']): continue
            if (kappa == s['kappa']):
                file_name = s['file_main'] + ".energies"
                f = open(file_name)
                data = (f.readlines()[-1].split())
                energy = float(data[1])
                f.close()
                eta = convert_to_number(s['eta'])
                if (eta_log): eta = math.log(eta)
                energy_graph[kappa].append((eta, energy))
        energy_graph[kappa].sort(key=get_key)

    files['kappas'] = []
                
    for kappakey in energy_graph:
        kappagraph = energy_graph[kappakey]
        #print(kappagraph)
        file_name = energy_dir + "kappa_" + kappakey + ".energies"
        f = open(file_name, 'w+')
        files['kappas'].append(file_name)
        for graph_point in kappagraph:
            eta = graph_point[0]
            energy = graph_point[1]
            f.write(str(eta))
            f.write(" ")
            f.write(str(energy))
            f.write('\n')
        f.close()
    
    return(files)



if (args.make_graphs): make_graphs()


def fill_active_sims(active_sims):
    flag = False
    i = 0
    while (i < float(args.max_active)):
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

def run_project():
    n_active_sims = 0
    i = 0
    active_sims = []
    while (i < float(args.max_active)):
        active_sims.append({})
        i += 1
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



if (args.run_project):
    run_project()
