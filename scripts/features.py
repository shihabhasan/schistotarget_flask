from Bio import SeqIO
from Bio.SeqUtils import ProtParam

import sys, os, subprocess
import re
import time, hashlib

from pepstat import pepstat
from garnier import garnier
from netcglyc import netcglyc
from netchop import netchop
from netoglyc import netoglyc
from netnglyc import netnglyc
from netphos import netphos
from tmhmm import tmhmm
from prop import prop
from signalp import signalp

def features(seqID, sequence):
    parameters=""
         
    seq=sequence.replace("X","").replace("*","").upper()

    X = ProtParam.ProteinAnalysis(seq)
    
    count=X.count_amino_acids()
    parameters=parameters+str(len(sequence))+","
    
    percent=X.get_amino_acids_percent()
    parameters=parameters+str(percent["A"])+","+str(percent["C"])+","+str(percent["D"])+","+str(percent["E"])+","+str(percent["F"])\
                   +","+str(percent["G"])+","+str(percent["H"])+","+str(percent["I"])+","+str(percent["K"])+","+str(percent["L"])+","+str(percent["M"])\
                   +","+str(percent["N"])+","+str(percent["P"])+","+str(percent["Q"])+","+str(percent["R"])+","+str(percent["S"])+","+str(percent["T"])\
                   +","+str(percent["V"])+","+str(percent["W"])+","+str(percent["Y"])+","

    

    parameters=parameters+str(X.molecular_weight())+","
    parameters=parameters+str(X.aromaticity())+","
    parameters=parameters+str(X.instability_index())+","
    parameters=parameters+str(X.isoelectric_point())+","
    parameters=parameters+str(X.gravy())+","
    parameters=parameters+str(X.secondary_structure_fraction()[0])+","+str(X.secondary_structure_fraction()[1])+","+str(X.secondary_structure_fraction()[2])+","
    parameters=parameters+str(X.molecular_weight()/len(seq))+","


    #-----------------------PEPSTATS-------------------------------
    
    inFileName=seqID.replace("|","_").replace(":","_")+"_"+hashlib.md5(time.asctime()).hexdigest()
    inFasta=open(inFileName, "w")
    inFasta.write(">"+seqID.replace("|","_")+"\n"+seq+"\n")
    inFasta.close()

    parameters=parameters+pepstat(inFileName)+","
                   

    tiny = count["A"]+count["C"]+count["G"]+count["S"]+count["T"]
    tiny_per = round(((float(tiny)/len(seq))*100),3)
    parameters=parameters+str(tiny_per)+","

    small = count["A"]+count["C"]+count["D"]+count["G"]+count["N"]+count["P"]+count["S"]+count["T"]+count["V"]
    small_per = round(((float(small)/len(seq))*100),3)
    parameters=parameters+str(small_per)+","

    aliphatic = count["I"]+count["L"]+count["V"]
    aliphatic_per = round(((float(aliphatic)/len(seq))*100),3)
    parameters=parameters+str(aliphatic_per)+","

    aromatic = count["F"]+count["H"]+count["W"]+count["Y"]
    aromatic_per = round(((float(aromatic)/len(seq))*100),3)
    parameters=parameters+str(aromatic_per)+","

    polar = count["D"]+count["E"]+count["H"]+count["K"]+count["N"]+count["Q"]+count["R"]+count["S"]+count["T"]
    polar_per = round(((float(polar)/len(seq))*100),3)
    parameters=parameters+str(polar_per)+","

    nonPolar = count["A"]+count["C"]+count["F"]+count["G"]+count["I"]+count["L"]+count["M"]+count["P"]+count["V"]+count["W"]+count["Y"]
    nonPolar_per = round(((float(nonPolar)/len(seq))*100),3)
    parameters=parameters+str(nonPolar_per)+","

    charged = count["D"]+count["E"]+count["H"]+count["K"]+count["R"]
    charged_per = round(((float(charged)/len(seq))*100),3)
    parameters=parameters+str(charged_per)+","
    
    acidic = count["D"]+count["E"]
    acidic_per = round(((float(acidic)/len(seq))*100),3)
    parameters=parameters+str(acidic_per)+","
    
    basic = count["H"]+count["K"]+count["R"]
    basic_per = round(((float(basic)/len(seq))*100),3)
    parameters=parameters+str(basic_per)+","
    


    #------------------------GARNIER--------------------------------

    parameters=parameters+garnier(inFileName)+","


    #--------------------------NetCGlyc--------------------------------

    parameters=parameters+netcglyc(inFileName)+","

    #--------------------------NetChop--------------------------------

    parameters=parameters+netchop(inFileName)+","
    

    #--------------------------NetNGlyc--------------------------------

    parameters=parameters+netnglyc(inFileName)+","
    

    #--------------------------NetPhos--------------------------------

    parameters=parameters+netphos(inFileName)+","


    #--------------------------TMHMM----------------------------------

    parameters=parameters+tmhmm(inFileName)+","


    #--------------------------ProP----------------------------------

    parameters=parameters+prop(inFileName)+","

    #--------------------------SignalP----------------------------------

    parameters=parameters+signalp(inFileName)

    os.remove(inFileName)
    
    return parameters

    
