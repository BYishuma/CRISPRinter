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
from Bio import SeqIO

from CRISPRinterLOCAL import CRISPRinterParams
from CRISPRinterLOCAL import OnlineResources
from CRISPRinterLOCAL import CRISPRinterResources

def main(spacerFile,repeatseq,proteinSeq,outDir,inunique,incov,inident,inrident,inmaxprotlen,blastdb,effectorflag):

    ### Now, reading protein length and check revise status
    proteinLength = len(proteinSeq)
    ProteinBlastDir = os.path.join(outDir,'ProteinBlast/')
    ProteinBlastResult = os.path.join(ProteinBlastDir,'ProteinBlastResult.txt')
    CRISPRinterParams.check_directory(ProteinBlastDir)

    ### Now, running protein blastp to predict related genomes
    ### Build temp protein fasta file
    tempProteinFile = os.path.join(ProteinBlastDir,'tempProteinSeq.fasta')
    CRISPRinterResources.build_temp_protein_fasta(proteinSeq,tempProteinFile)
    CRISPRinterResources.run_local_blastp(tempProteinFile,ProteinBlastResult,blastdb,effectorflag)
    os.system("rm "+tempProteinFile)
    if CRISPRinterParams.check_null_files(ProteinBlastResult):
        print("No related protein found. Spacer revise failed.")
        return False
        
    ### Get significant hits and return putative source tax name
    SignificantHitResult = os.path.join(ProteinBlastDir,"SignificantHits.txt")
    LevelBox = CRISPRinterResources.get_local_significant_blastp_hits(ProteinBlastResult,SignificantHitResult,incov,inident)
    if CRISPRinterParams.check_null_files(SignificantHitResult):
        print("No significant hits found. Spacer revise failed. Use original spacer file")
        return False

    ### Now, select the correct DRs based on the MinCED result and repeat file
    CasFullInfoPickle = os.path.join(blastdb,"CasInfoPICKLE/",effectorflag+".pickle")
    spacerPickle = os.path.join(blastdb,"SpacerCollection.pkl")
    SelectedSpacerFilesDict = CRISPRinterResources.select_correct_dr(spacerPickle,repeatseq,inrident,CasFullInfoPickle,LevelBox)

    ### Now, check the selected spacer file dict, if the dict is empty, return False
    if len(SelectedSpacerFilesDict) == 0:
        print("No correct DR found. Spacer revise failed. Use original spacer file")
        return False
    
    ### Now, revise the spacer file
    RevisedSpacerFile = os.path.join(outDir,'RevisedSpacer.txt')
    #CRISPRinterResources.revise_spacer_file(spacerFile,MinCEDDir,SelectedSpacerFilesDict,RevisedSpacerFile,inunique)
    CRISPRinterResources.revise_spacer_file(spacerFile,SelectedSpacerFilesDict,inunique,RevisedSpacerFile)

    ### Return revised spacer file
    return RevisedSpacerFile

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])