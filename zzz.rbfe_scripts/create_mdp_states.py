import argparse

parser =  argparse.ArgumentParser()
parser.add_argument("-f", "--mdpfile", help="Gromacs free energy MDP file to be modified",
                    type=str, nargs=1, required=True)

args = parser.parse_args()
if len(args.mdpfile) != 1:
    parser.print_help()

string1 = "fep-lambdas ="
string2 = "init-lambda-state ="
mdpfile = args.mdpfile[0]

with open(mdpfile, "r+") as mdp:
    lines = mdp.readlines()

for line in lines:
    if string1 in line:
        tmp = line.split(" ")
        # remove the first two elements
        tmp.pop(0)
        tmp.pop(0) 
        num_lam = len(tmp)
    else: continue
print(f"There are {num_lam} states")

for x in range(0, num_lam):
    tmp = mdpfile.split(".X.")
    new_mdpfile = tmp[0] + f".{x}.mdp"
    g = open(new_mdpfile, "w+")
    for line in lines:
        if string2 not in line:
            g.write(line)
        else:
            newline = f"init-lambda-state = {x} \n"
            g.write(newline)
    g.close()

print(f"All states' MDP files for {mdpfile} are ready.")

