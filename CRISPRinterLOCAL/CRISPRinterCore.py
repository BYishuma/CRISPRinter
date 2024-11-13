#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
CRISPRinterLOCAL - A CRISPR homologous enhancement-based host/phage interaction prediction tool, LOCAL
Author: 
BNU-HKBU UIC
Email: 
'''

import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import os
import time

print("===============================================================")
print("Start loading resources...")
resourcesTime = time.time()
runningTime = time.time()

from CRISPRinterLOCAL import CRISPRinterParams
from CRISPRinterLOCAL import CRISPRinterReviser
from CRISPRinterLOCAL import CRISPRinterResources
from CRISPRinterLOCAL import OnlineResources

print("Loading resources finished. Time used: %s seconds" % (time.time() - resourcesTime))

def main():
    print("===============================================================")
    print("Initializing parameters...")

    arg_parser = CRISPRinterParams.get_arg_parser()
    opts = arg_parser.parse_args()

    InputGenome = opts.input
    OutDir = opts.outdir
    BlastMode = opts.blastmode
    BacteriaDatabase = opts.bacteriaDB
    PhageDatabase = opts.phageDB
    ProteinCOV = opts.pcovs
    ProteinIDENT = opts.pident
    UniqueMode = opts.unique
    RepeatIDENT = opts.rident
    MaxProteinNum = opts.MaxProteinNum
    CRISPRCasTyperDB = opts.cctyperDB

    ### FILE CHECKS

    CRISPRinterParams.check_environment()
    CRISPRinterParams.check_outdir(OutDir)

    ### MODE CHECKS AND SPACER REVISE

    print("Initializing parameters finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    ### Predict all potential CRISPR-Cas systems of the input bacteria genome
    CRISPRCasTyperDir = os.path.join(OutDir,"CRISPRCasTyperResult/")
    CRISPRinterResources.run_CRISPRCasTyper(InputGenome,CRISPRCasTyperDB,CRISPRCasTyperDir)

    CRISPRinterResources.check_CRISPRCas_status(CRISPRCasTyperDir)
    CRISPRinfoLabeledDict = CRISPRinterResources.get_CRISPRCas_info(CRISPRCasTyperDir)

    ### [crisprid]:[ContigID,SpacerFile,RepeatSeq,ProteinSeq,Label,CasFlag]
    ### Label: WithProt & NoProt
    ### CasFlag: Cas9/Cas3/Cas12/Cas13/None ...

    ReviserDir = os.path.join(OutDir,"Reviser/")
    CRISPRinterParams.check_directory(ReviserDir)

    for CRISPRunit, CRISPRunitInfo in CRISPRinfoLabeledDict.items():
        CRISPRLabel = CRISPRunitInfo[4]
        CasFlag = CRISPRunitInfo[5]
        # print("This is CasFlag: "+CasFlag)
        targetDir = os.path.join(ReviserDir,CRISPRunit)
        if CRISPRLabel == "WithProt":
            ### Start CRISPR homologous enhancement

            CRISPRinterParams.check_directory(targetDir)
            finalspacer = CRISPRinterReviser.main(CRISPRunitInfo[1],CRISPRunitInfo[2],CRISPRunitInfo[3],targetDir,UniqueMode,ProteinCOV,ProteinIDENT,RepeatIDENT,MaxProteinNum,BacteriaDatabase,CasFlag)
            if finalspacer == False:
                finalspacer = CRISPRunitInfo[1]
        else:
            ### No enhancement, use spacer file to predict phage
            finalspacer = CRISPRunitInfo[1]
        ### Now, modify the finalspacer file

        spaceridLenDict = CRISPRinterResources.modify_spacer_file(finalspacer)
        ### Now, run spacer blast
        ### SETTING SPACER BLAST DIRECTORY
        SpacerBlastDir = os.path.join(targetDir,"SpacerBlast/")
        CRISPRinterParams.check_directory(SpacerBlastDir)
        SpacerInputDir = os.path.join(SpacerBlastDir,"Input/")
        CRISPRinterParams.check_directory(SpacerInputDir)
        SpacerBlastOutputDir = os.path.join(SpacerBlastDir,"RawOutput/")
        CRISPRinterParams.check_directory(SpacerBlastOutputDir)

        ### CHECK SPACER LENGTH, IF LENGTH > 20, CUT TO 20
        SpacerNumbers = os.popen("grep -c '>' %s" % finalspacer).read().strip()
        if int(SpacerNumbers) > 20:
            CRISPRinterResources.split_spacer_file(finalspacer,SpacerInputDir,int(SpacerNumbers))
        else:
            os.system("cp %s %s" % (finalspacer,SpacerInputDir))

        ### BLAST SPACER AGAINST REPEAT DATABASE
        for InputSpacerFile in os.listdir(SpacerInputDir):
            InputSpacerFile = os.path.join(SpacerInputDir,InputSpacerFile)
            OutputSpacerFile = os.path.join(SpacerBlastOutputDir,os.path.basename(InputSpacerFile))
            OnlineResources.run_spacer_blast(InputSpacerFile,OutputSpacerFile,BlastMode)

        ### MERGE BLAST OUTPUT IF THERE HAVE MULTIPLE BLAST OUTPUT
        MergedSpacerBlastOutput = os.path.join(SpacerBlastDir,"MergedSpacerBlastOutput.txt")
        if len(os.listdir(SpacerBlastOutputDir)) > 1:
            CRISPRinterResources.merge_blast_output(SpacerBlastOutputDir,MergedSpacerBlastOutput)
        else:
            os.system("cp %s %s" % (os.path.join(SpacerBlastOutputDir,os.listdir(SpacerBlastOutputDir)[0]),MergedSpacerBlastOutput))

        ### GET SIGNIFICANT SPACER BLAST OUTPUT AND NON-SIGNIFICANT HIT SPACER ID
        SignificantSpacerBlastOutput = os.path.join(SpacerBlastDir,"SignificantSpacerBlastOutputRaw.txt")
        InSignificantSpacerID = CRISPRinterResources.get_significant_spacer_blast_output(MergedSpacerBlastOutput,SignificantSpacerBlastOutput,spaceridLenDict)
        
    ### SET FINAL SPACER HIT FILE
    FinalSpacerHitFile = os.path.join(OutDir, "FinalSpacerHit.txt")
    CRISPRinterResources.merge_final_spacer_hit(ReviserDir, FinalSpacerHitFile)
    FinalSpacerHitTaxNameFile = os.path.join(OutDir, "FinalSpacerHitTaxNameFile.txt")
    OnlineResources.subject_id_2_tax_name(FinalSpacerHitFile, FinalSpacerHitTaxNameFile)
    Result1 = os.path.join(OutDir, "taxid_count.txt")
    Result2 = os.path.join(OutDir, "Speices.txt")
    OnlineResources.final_prediction(FinalSpacerHitTaxNameFile, Result1, Result2)
        
if __name__ == "__main__":
    main()
