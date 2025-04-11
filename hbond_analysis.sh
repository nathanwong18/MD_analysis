#!/bin/bash
#Automate hbond analysis over a set of trajectories, start with just one 
#
#Given a .gro file, this script will:
#
#-Identify a 6 angstrom cutoff region around a charge carrying cofactor
#-Identify all residues within that region
#-Identify hbonding pair distances within 3.5 angstroms and 30 degrees (be able to change this parameter)
#-Could use GMX or VMD hbonding modules need to iterate over all frames...why not both and double check for redundency
#-output to txt/dat file
#
#1 microsecond of simulation
#
#
#Amino acid library
#Positive AA
POS=(ARG,LYS,HIS) #Positive AA
 #Arginine (ARG)
 #Lysine (LYS)
 #Histidine (HIS)
#Negative AA
NEG=(GLU,ASP) #Negative AA
 #Glutamic acid (GLU)
 #Aspartic acid (ASP)
#Polar AA
POL=(SER,THR,ASN,GLN,CYS) #Polar AA
 #Serine (SER)
 #Threonine (THR)
 #Asparagine (ASN)
 #Glutamine (GLN)
 #Cysteine (CYS)
#Non-polar AA
NP=(ALA,VAL,LEU,ILE,PRO,MET,PHE,TRP,TYR) #Non-polar AA
 #Alanine (ALA)
 #Valine (VAL)
 #Leucine (LEU)
 #Isoleucine (ILE)
 #Proline (PRO)
 #Methionine (MET)
 #Phenylalanine (PHE)
 #Tryptophan (TRP)
 #Tyrosine (TYR) 
#
#Load the VMD module
echo "module load vmd/1.9.3-intel" > hbond.txt
#
#Parameters/Inputs
cutoff=3.5
angle=30
gro=md_01.gro #INSERT YOUR GRO FILE HERE
pdb=oo-chains.pdb #INSERT YOUR PDB FILE HERE
xtc=small.xtc #INSERT YOUR XTC FILE HERE
energy_file=vertical_00-OO_01-OR.dat
cofactor=(500) # 500 resid of [2Fe-2S] cluster index residue 68, 144
end=10


echo "Starting hbond analysis using VMD"
#VMD
for i in "${cofactor[@]}"; do
    cofactorID=$i
    #1.Generate a txt file to run on VMD
        echo "package require hbonds" >> hbond.txt
        echo "mol load pdb $pdb" >> hbond.txt
        echo "mol addfile $xtc type xtc waitfor all" >> hbond.txt
        echo "set cutoff $cutoff" >> hbond.txt
        echo "set angle $angle" >> hbond.txt
        #echo "set traj [mol new .gro]" >> hbond.txt
        echo "set num_steps [molinfo top get numframes]" >> hbond.txt
        #echo "set end [expr $num_steps-1]" >> hbond.txt
        echo "set end $end" >> hbond.txt
        echo "set sel1 [atomselect top \"not water and same residue as within $cutoff of (resid $cofactorID and chain 1)\"]" >> hbond.txt
        echo "set sel2 [atomselect top \"not water and same residue as within $cutoff of (resid $cofactorID and chain 1)\"]" >> hbond.txt
        #echo "hbonds -sel1 \$sel1 -sel2 \$sel2 -frames 0:\$end -dist 3.5 -ang 60 -polar yes -plot no -writefile yes -outfile ${cofactorID}_Hbond.dat  -type unique -detailout ${cofactorID}_Hbond-detailed.dat" >> hbond.txt
        echo "hbonds -sel1 \$sel1 -sel2 \$sel2 -frames 0:\$end -dist 3.5 -ang 60 -polar yes -plot no -writefile yes -outfile ${cofactorID}_Hbond.dat  -type unique -detailout ${cofactorID}_Hbond-detailed.dat" >> hbond.txt
        echo "quit" >> hbond.txt
    #2.Run VMD
        vmd -dispdev text -e hbond.txt
done 
echo "hbond analysis complete"
#3.Calculate distance of H-bonding pairs identified from previous step from VMD
echo "processing hbond-detailed.dat"
#process the hbond-detailed.dat file
#This file contains the following columns:
head -n 1 ${cofactorID}_Hbond-detailed.dat
awk 'BEGIN {
    # Print header row
    printf "%-12s %-13s %-14s %-13s %-12s %-13s %-14s %-13s\n", "Donor_AA", "Donor_Num", "DonorType", "DonorAtom", "Acceptor_AA", "Acceptor_Num", "AcceptorType", "AcceptorAtom"
    print "-------------------------------------------------------------------------------------------------------------------"
}
NR > 2 {
    # Split donor and acceptor fields using dash and whitespace as delimiters
    nd = split($1, d, /[- ]+/)
    na = split($2, a, /[- ]+/)
    
    # For donor field: extract three-letter code and number
    donorAA = substr(d[1], 1, 3)
    donorNum = ""
    if(match(d[1], /[0-9]+/))
        donorNum = substr(d[1], RSTART, RLENGTH)
    
    # Concatenate the remaining donor parts (if any)
    donorType = (d[2] ? d[2] : "")
    donorAtom = (d[3] ? d[3] : "")

    # For acceptor field: extract three-letter code and number
    acceptorAA = substr(a[1], 1, 3)
    acceptorNum = ""
    if(match(a[1], /[0-9]+/))
        acceptorNum = substr(a[1], RSTART, RLENGTH)
    
    # Concatenate the remaining acceptor parts (if any)
    acceptorType = (a[2] ? a[2] : "")
    acceptorAtom = (a[3] ? a[3] : "")
    # Print in formatted columns
    printf "%-12s %-13s %-14s %-13s %-12s %-13s %-14s %-13s\n", donorAA, donorNum, donorType, donorAtom, acceptorAA, acceptorNum, acceptorType, acceptorAtom
}' ${cofactorID}_Hbond-detailed.dat > hbondnetwork.dat
echo "hbond-detailed.dat processed"
#This file contains the following columns:
#HDonor: Hydrogen bond donor
#HAcceptor: Hydrogen bond acceptor
#HbondDistance: Distance between donor and acceptor
#HbondAngle: Angle between donor, hydrogen, and acceptor
#RedoxPotential: Redox potential of the donor-acceptor pair
#VE1: Redox state of ground state (ie OO)
#VE2: Redox state of excited state (ie OR, RO, RR)
#Table Format
#     HDonor  HAcceptor HbondDistance HbondAngle RedoxPotential   VEgs   VEes   Mutation  Type     WithinHbondcutoff (Y/N 0/1)" #create header
#ie    ARG73   FMR500        3.5         30           0.0          OO     RO      WT      Main         0

#2. Read hbondnetwork.dat for Hbond pairs and calculate distances for each frame 
#       -start with the first pair then iterate over the rest
#       -Set VMD atomselect to the donor and acceptor atoms and resid
echo "Reading hbondnetwork.dat"
j=0

while IFS= read -r line; do
    # Skip the first two lines
    if [[ $((++line_num)) -le 2 ]]; then
        continue
    fi
    echo "Processing line $line_num: $line"
    # Read the donor and acceptor from the line
    donorAtom=$(echo "$line" | awk '{print $4}')
    acceptorAtom=$(echo "$line" | awk '{print $8}')
    # Read the donor and acceptor numbers from the line
    donorNum=$(echo "$line" | awk '{print $2}')
    acceptorNum=$(echo "$line" | awk '{print $6}')
    echo "Donor Atom: $donorAtom"
    echo "Acceptor Atom: $acceptorAtom"
    echo "Donor Number: $donorNum"
    echo "Acceptor Number: $acceptorNum"
    echo "nextline"
    touch VMDHbondout.txt
    VMDout=VMDHbondout.txt
    echo "mol load pdb $pdb" > hbond_distance.txt
    echo "mol addfile $xtc type xtc waitfor all" >> hbond_distance.txt
    echo "set D [[atomselect top \"name $donorAtom and resid $donorNum and chain 1\"] get index]" >> hbond_distance.txt
    echo "set A [[atomselect top \"name $acceptorAtom and resid $acceptorNum and chain 1\"] get index]" >> hbond_distance.txt
    echo "set d1 \"\$D \$A\"" >> hbond_distance.txt
    echo "set nf [molinfo top get numframes]" >> hbond_distance.txt
    echo "set outfile [open $VMDout w]" >> hbond_distance.txt
    echo "set DA [measure bond \$d1 frame all]" >> hbond_distance.txt
    echo "set i 0" >> hbond_distance.txt
    echo "foreach x \$DA {" >> hbond_distance.txt
    echo "    incr i" >> hbond_distance.txt
    #echo "    puts \"frame \$i of \$nf\" " >> hbond_distance.txt
    echo "    puts \$outfile \"\$i \$x\" " >> hbond_distance.txt
    echo "}" >> hbond_distance.txt
    echo "close \$outfile" >> hbond_distance.txt
    echo "quit" >> hbond_distance.txt
    # Run VMD to calculate the distance
    vmd -dispdev text -e hbond_distance.txt < /dev/null
    ## Read the distance from the output
    #NOTE distance is from donor to acceptor average NH bond distance is 1.01 A
    #NOTE distance is from acceptor to donor average OH bond distance is 0.95 A
    # Append the distance to the hbondnetwork.dat file
    echo "Distance calculated"
#3. Combine all compiled data into a single file
    if [ $j -eq 0 ]; then
        touch output.txt
        header1=$(head -n 1 hbondnetwork.dat) 
        echo "$header1 "       Hbond Distance  "        "    LJ  "     "   Coulomb "  " "Potential" " " > output.txt
        header2=$(sed -n '2p' hbondnetwork.dat)
        echo "$header2 ----------------" >> output.txt
        j=1
        echo "Header added to output.txt"
    fi 
    if [ -f VMDHbondout.txt ]; then
    echo "Appending data to output.txt"
    line_count=$(wc -l < VMDHbondout.txt)
    echo "FRAME INFO: VMDHbondout.txt has $line_count lines."
    echo "Save the first $line_count lines in $energy_file"
    paste VMDHbondout.txt <(head -n "$line_count" $energy_file) > temp.txt  
    echo "Adding Distance and energy"
    while IFS= read -r m; do
        distanceAndEnergy=$(echo "$m" | awk '{ printf "%-26s %-9s %-9s %-9s", $2, $5, $6, $7}')
        echo "$line $distanceAndEnergy" >> output.txt
    done < temp.txt
    rm temp.txt
    else
        echo "Error: VMDHbondout.txt not found. Skipping distance calculation."
    fi
echo "\$------------------------------------------------------------------------------------------------"
echo ""
echo ""
done < hbondnetwork.dat

