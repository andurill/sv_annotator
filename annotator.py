'''
Copyright (c) 2019 Memorial Sloan-Kettering Cancer Center.

This script is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 2.1 of the License, or
any later version.

This script is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
documentation provided hereunder is on an "as is" basis, and
Memorial Sloan-Kettering Cancer Center has no obligations to provide
maintenance, support, updates, enhancements or modifications.  In no
event shall Memorial Sloan-Kettering Cancer Center be liable to any
party for direct, indirect, special, incidental or consequential damages,
including lost profits, arising out of the use of this software and its
documentation, even if Memorial Sloan-Kettering Cancer Center has been
advised of the possibility of such damage.  See the GNU Lesser General
Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.

Created on Jan 20, 2019

@author: balakra1@mskcc.org,jayakumg@mskcc.org
@version: 1.0.0

'''

from bp import *
from constants.constants import *
from utils import *
import re

class SVAnnotator(object):
    """
    SV annotation class

    """

    def __init__(self,gene1,gene2,chr1,chr2,pos1,pos2,desc1,desc2,svType,event,logger):
        r"""
        Constrcutor for SV annotation class

        Parameters
        ----------
        self
        gene1 : string
        gene1 : string
        chr1 : int
        chr2 : int
        pos1 : int
        pos2 : int
        desc1 : string
        desc2 : string
        svType : string
        event : string
        logger : logger

        Returns
        -------
        None

        Raises
        -------
        Exception
            Critical exception while getting annotation parameters

        Notes
        -------
        calls __annotator fuction
        """
        self.__logger = logger
        self.__gene_1 = gene1
        self.__chr_1 = chr1
        self.__pos_1 = pos1
        self.__desc_1 = desc1
        self.__gene_2 = gene2
        self.__chr_2 = chr2
        self.__pos_2 = pos2
        self.__desc_2 = desc2
        self.__is_oncogenic = ""
        self.__sv_type = SV_TYPE_EXP.get(svType,"UNKNOWN")
        self.__sv_type_code = svType
        self.__event_info = event
        self.__sv_special_type = SV_SPECICAL_TYPE_UNK
        self.__is_known_oncokb_fusion = False
        self.__initialize_notes()
        self.__annotator()
        # try:
        #     self.__annotator()
        # except Exception as e:
        #     self.__logger.error(e)
        #     self.__logger.error("Exception while trying to annotate SV")
        #     self.__logger.error(gene1,gene2,chr1,chr2,pos1,pos2,desc1,desc2,svType,event)
        #     raise Exception(e)

    def __initialize_notes(self):
        r"""
        Initialize Notes
        """
        self.__annotation_note = self.__comment_note = self.__bkpsites = self.__gene_note = self.__sig_note = \
        self.__fusion_type = self.__transcript_coordinates = self.__transcript_coordinates_postfix = self.__annotation_note_sv_type = \
        self.__comment_gene_note = self.__note_prefix = self.__misc_note = self.__exons_note = ""

    def __shuffle(self,__t1,__t2):
        r"""
        Shuffle two objects

        Parameters
        ----------
        self
        __t1 : object1
        __t2 : object2

        Returns
        -------
        object2, object1
        """
        return __t2,__t1

    def __join(self,__l):
        r"""
        Concatnate list of string

        Parameters
        ----------
        self
        __l : list

        Returns
        -------
        string
        """
        return ' '.join(filter(None, __l))

    def __set_is_oncokb_fusion(self):
        r"""
        Check oncoKb API to see if fusion of genes is Oncgenic

        Parameters
        ----------
        self

        Returns
        -------
        None

        Other Parameters
        ----------------
        __logger : logger adapter

        Notes
        -----
        This function call the is_oncokb_fusion function in utils and 
        sets __is_known_oncokb_fusion based on return value

        """
        self.__is_known_oncokb_fusion = is_oncokb_fusion(self,self.__logger)

    def __set_fusion_partners(self):
        r"""
        Set SV gene order based on gene order in `__event_info`

        Parameters
        ----------
        self

        Returns
        -------
        None

        Raises
        ----------------
        Exception
            Genes format in __event_info does not match format {gene1:gene2}
        Exception
            Genes in __event_info does not match __gene_1 & __gene_2
        """
        try:
            __fusion_gene_1 , __fusion_gene_2 = self.__event_info[self.__event_info.find("{")+1:self.__event_info.find("}")].split(":")
        except Exception as e:
            self.__logger.error("Incorrect fusion gene format")
            self.__logger.error(e)
            raise Exception("Incorrect fusion gene format")

        if __fusion_gene_2 == self.__partner_1.get_gene() and __fusion_gene_1 == self.__partner_2.get_gene():
            self.__partner_1, self.__partner_2 = self.__shuffle(self.__partner_1, self.__partner_2)
        elif __fusion_gene_1 == self.__partner_1.get_gene() and __fusion_gene_2 == self.__partner_2.get_gene(): pass
        else:
            self.__logger.error("Fusion Partner and Genes do not match")
            raise Exception("Fusion Partner and Genes do not match")

    def __set_annotation_gene_note(self,__bkp_1,__bkp_2=None):   
        r"""
        Set gene_note part of Annotation & Comment

        Parameters
        ----------
        self
        __bkp_1 : Breakpoint
            first breakpoint or partner
        __bkp_2 : Breakpoint
            second breakpoint or partner
            default `None`

        Returns
        -------
        None

        Notes
        ----------------
        sets `self.__gene_note` & `self.__comment_gene_note`

        """
        if __bkp_2 != None and __bkp_1.get_gene() != __bkp_2.get_gene():
            self.__gene_note = "{bk1} - {bk2}".format(bk1=self.__set_annotation_gene_note_for_single_gene(__bkp_1),bk2=self.__set_annotation_gene_note_for_single_gene(__bkp_2))
            self.__comment_gene_note = "{g1} - {g2}".format(g1=__bkp_1.get_gene(),g2=__bkp_2.get_gene())
        else:
            self.__gene_note = self.__set_annotation_gene_note_for_single_gene(__bkp_1)
            self.__comment_gene_note = __bkp_1.get_gene()

    def __set_annotation_gene_note_for_single_gene(self,__bkp):
        r"""
        Form gene_note part of Annotation for a single gene

        Parameters
        ----------
        self
        __bkp : Breakpoint
            breakpoint or partner

        Returns
        -------
        `gene_note` : format <gene> (<Refseq Transcript>)

        Notes
        ----------------
        Refseq Transcript is set in the Breakpoint class.

        """
        return "{gene} ({tx})".format(gene=__bkp.get_gene(),tx=__bkp.get_transcript())

    def __set_annotation_tx_coord_note(self,__bkp_1,__bkp_2=None):
        r"""
        Set transcript_coordinates part of the notes

        Parameters
        ----------
        self
        __bkp_1 : Breakpoint
            first breakpoint or partner
        __bkp_2 : Breakpoint
            second breakpoint or partner
            default `None`

        Returns
        -------
        None

        Notes
        ----------------
        sets `self.__transcript_coordinates`
        """
        if __bkp_2 != None:
            if __bkp_1.get_gene() != __bkp_2.get_gene():
                self.__transcript_coordinates = "{bk1}_{bk2}".format(bk1=self.__set_annotation_tx_coord_note_for_single_gene(__bkp_1),bk2=self.__set_annotation_tx_coord_note_for_single_gene(__bkp_2))
            else:
                self.__transcript_coordinates = "{cdna1}_{cdna2}".format(cdna1=__bkp_1.get_cdna(),cdna2=__bkp_2.get_cdna())
        else:
            self.__transcript_coordinates = self.__set_annotation_tx_coord_note_for_single_gene(__bkp_1)

    def __set_annotation_tx_coord_note_for_single_gene(self,__bkp):
        r"""
        Form transcript_coordinates part of the notes for a single gene

        Parameters
        ----------
        self
        __bkp : Breakpoint
            breakpoint or partner

        Returns
        -------
        transcript_coordinates_for_single_gene

        Notes
        ----------------
        format <cdna>:<gene> if breakpoint is in coding else <cdna>

        """
        if __bkp.get_is_coding():
            return "{cdna}:{gene}".format(cdna=__bkp.get_cdna(),gene=__bkp.get_gene())
        else:
            return "{cdna}".format(cdna=__bkp.get_cdna())

    def __set_annotation_tx_coord_note_postfix(self,__bkp):
        r"""
        Sets __transcript_coordinates_postfix which is tail end of annotation

        Parameters
        ----------
        self
        __bkp : Breakpoint
            breakpoint or partner

        Returns
        -------
        None

        Notes
        ----------------
        sets `__transcript_coordinates_postfix` 
        format `_chr<chromosome>:g.<position>`

        """
        self.__transcript_coordinates_postfix = "_chr{chr}:g.{pos}".format(chr=__bkp.get_chrom(),pos=__bkp.get_pos())

    def __set_annotation_tx_coord_note_for_tra(self):
        r"""
        Set transcript_coordinates part of the a Translocation SV

        Parameters
        ----------
        self

        Returns
        -------
        None

        Notes
        ----------------
        sets `self.__transcript_coordinates` for Translocation SV based on chromosome order
        """
        if self.check_if_sv_type_translocation():
            if (self.__partner_1.get_chrom().isdigit() and self.__partner_2.get_chrom().isdigit() and int(self.__partner_2.get_chrom()) < int(self.__partner_1.get_chrom())) or \
                (self.__partner_1.get_chrom() in ("X","Y") and self.__partner_2.get_chrom().isdigit() ) or (self.__partner_1.get_chrom() == "Y" and self.__partner_1.get_chrom() == "X"):
                self.__transcript_coordinates = get_coordinate(self.__partner_2,self.__partner_1,self.__logger)
            else: self.__transcript_coordinates = get_coordinate(self.__partner_1,self.__partner_2,self.__logger)

    def __generate_note_prefix(self):
        r"""
        Forms the prefix part of the comments based on SV type

        Parameters
        ----------
        self

        Returns
        -------
        __prefix : string
            Prefix string for comment

        """
        __prefix = ''
        try:
            __conj =  NOTE_CONNECTION_VOWELS if self.check_if_sv_type_inversion() or self.check_is_intragenic() else NOTE_CONNECTION_NON_VOWELS
            __sv_type = self.__sv_type.lower()
            if self.__is_known_oncokb_fusion:
                __prefix = NOTE_CONNECTION_ONCOKB_FUSION
            elif self.get_is_protein_fusion():
                __prefix = "{conn} {sv_type} {note}".format(conn=__conj, sv_type=self.__sv_type.lower(), note=NOTE_CONNECTION_PROTEIN_FUSION)
            elif self.check_is_intragenic():
                if self.__partner_1.get_intron() == self.__partner_2.get_intron() != None:
                    __prefix = "{conn} {intra} {sv_type} {bp}".format(conn=__conj, intra=NOTE_CONNECTION_INTRAGENIC, sv_type=self.__sv_type.lower(), bp=NOTE_CONNECTION_TRA_BOTH_BP)
                else:
                    __prefix = "{conn} {intra} {sv_type} {of}".format(conn=__conj, intra=NOTE_CONNECTION_INTRAGENIC, sv_type=self.__sv_type.lower(), of=NOTE_CONNECTION_OF)
            elif self.check_if_sv_type_translocation():
                __bp = NOTE_CONNECTION_TRA_ONE_BP
                if all([self.__partner_1.get_is_panel_and_coding(), self.__partner_2.get_is_panel_and_coding()]):
                    __bp = NOTE_CONNECTION_TRA_BOTH_BP
                __prefix = "{conn} {sv_type} {bp}".format(conn=__conj, sv_type=self.__sv_type.lower(), bp=__bp)
            else:
                __prefix = "{conn} {sv_type} {of}".format(conn=__conj, sv_type=self.__sv_type.lower(), of=NOTE_CONNECTION_OF)
        except Exception as e:
            self.__logger.error("Exception in __generate_note_prefix")
            self.__logger.error(e)
        finally:
            return __prefix

    def overwrite_protein_fusion(self):
        r"""
        This function overrides the oncokb sv type based on hotspot, kinase and entire_kinase values

        Parameters
        ----------
        self

        Returns
        -------
        __prefix : string
            Prefix string for comment

        Notes
        ----------------
        Sets __fusion_type to `fusion` for known oncokb fusion or protein fusion where any of partners is
        hotspot and entire kinase is invloved. 
        Overwrites __fusion_type to `rearrangement` if any of partners is kinase and entire kinase is 
        not invloved

        """
        if not self.__is_known_oncokb_fusion :
            if (not self.__partner_1.get_is_hotspot() and self.__partner_1.get_is_kinase()) or \
                (not self.__partner_2.get_is_hotspot() and self.__partner_2.get_is_kinase()):
                pass
            elif self.get_is_protein_fusion() and \
                any([all([self.__partner_1.get_is_hotspot(),self.__partner_1.check_is_entire_kinase_involved()]),all([self.__partner_2.get_is_hotspot(),self.__partner_2.check_is_entire_kinase_involved()])]):
                self.__is_known_oncokb_fusion = True 
        elif any([all([self.__partner_1.get_is_kinase_hotspot(),not self.__partner_1.check_is_entire_kinase_involved()]),all([self.__partner_2.get_is_kinase_hotspot(),not self.__partner_2.check_is_entire_kinase_involved()])]):
            self.__is_known_oncokb_fusion = False 
        self.__fusion_type = NOTE_FUSION_TYPE_FUSION if self.__is_known_oncokb_fusion else NOTE_FUSION_TYPE_REARRANGEMENT


    def __check_fusion_genes_in_coding(self,__bkp_1,__bkp_2):
        r"""
        This function checks if the fusion genes is in coding and sets transcript
        coordinate note postfix

        Parameters
        ----------
        self
        __bkp_1 : Breakpoint
            first breakpoint or partner
        __bkp_2 : Breakpoint
            second breakpoint or partner
            default `None`

        Returns
        -------
        __bkp_1 : Breakpoint
        __bkp_2 : Breakpoint

        Notes
        ----------------
        Sets the __set_annotation_tx_coord_note_postfix based on breakpoints
        If one of the breakpoints is not in coding its set to None and coding breakpoint
        is set as first breakpoint

        """
        try:
            if not __bkp_2.get_is_coding():
                self.__set_annotation_tx_coord_note_postfix(__bkp_2)
                __bkp_2 = None
            if not __bkp_1.get_is_coding():
                self.__set_annotation_tx_coord_note_postfix(__bkp_1)
                __bkp_1 = __bkp_2
                __bkp_2 = None
        except Exception as e:
            self.__logger.error("Execption in __check_fusion_genes_in_coding")
            self.__logger.error(e)
        finally:
            return __bkp_1,__bkp_2

    def __generate_comment_notes(self,__bkp_1,__bkp_2):
        r"""
        This function generates exon_note, misc_note, sig_note, note_prefix of comments

        Parameters
        ----------
        self
        __bkp_1 : Breakpoint
            first breakpoint or partner
        __bkp_2 : Breakpoint
            second breakpoint or partner
            default `None`

        Returns
        -------
        None

        """
        try:
            self.__exons_note = get_exons_involved(self,__bkp_1,__bkp_2,self.__logger)
            self.__misc_note = get_misc_notes(self)
            self.__exons_note = get_exons_involved(self,__bkp_1,__bkp_2,self.__logger)
            self.__sig_note = (functional_significance(self) or "")
            self.__note_prefix = self.__generate_note_prefix()
        except Exception as e:
            self.__logger.error("Execption while forming notes")
            self.__logger.error(e)

    def __generate_note_variables(self,__bkp_1,__bkp_2):
        r"""
        This is the main function that calls various functions to get partners and  partener order 
        and sets different parts of the annotation & comment

        Parameters
        ----------
        self
        __bkp_1 : Breakpoint
            first breakpoint or partner
        __bkp_2 : Breakpoint
            second breakpoint or partner
            default `None`

        Returns
        -------
        None

        """
        try:
            self.__fusion_type = NOTE_FUSION_TYPE_REARRANGEMENT
            self.__set_is_oncokb_fusion()

            if self.get_is_protein_fusion(): 
                self.__set_fusion_partners()
                __bkp_1,__bkp_2 =  self.__check_fusion_genes_in_coding(self.__partner_1,self.__partner_2)
            elif self.__check_if_both_breakpoints_in_panel_coding():
                if self.check_is_intragenic() and self.__partner_1.get_strand() == NEGATIVE_STRAND_IDENTIFIER: 
                    __bkp_1,__bkp_2 = self.__partner_1, self.__partner_2 = self.__shuffle(self.__partner_1, self.__partner_2)
            elif self.__partner_1.get_is_panel_and_coding():
                self.__set_annotation_tx_coord_note_postfix(__bkp_2)
                __bkp_2 = None
            else:
                self.__set_annotation_tx_coord_note_postfix(__bkp_1)
                __bkp_1,__bkp_2 = self.__partner_2, None

            self.__set_annotation_gene_note(__bkp_1,__bkp_2)
            self.__set_annotation_tx_coord_note(__bkp_1,__bkp_2)
            self.__generate_comment_notes(__bkp_1,__bkp_2)

            if self.check_if_sv_type_translocation():
                self.__set_annotation_tx_coord_note_for_tra()
                self.__transcript_coordinates_postfix = ""
            else:
                self.__annotation_note_sv_type = self.__sv_type_code.lower()

        except Exception as e:
            self.__logger.error("Execption while forming annotation")
            self.__logger.error(e)

    def __annotator(self):
        r"""
        This is the primary functiom which forms the breakpoint instance from input arguments
        and generates the variable component for notes

        Parameters
        ----------
        self

        Returns
        -------
        None

        Raises
        ----------------
        Exception
            Genes not in panel
        Exception
            Both breakpoints not in coding

        """
        try:
            __bkp_1 = self.__partner_1 = Breakpoint(self.__chr_1, self.__pos_1, self.__gene_1, self.__desc_1, self.__logger)
            __bkp_2 = self.__partner_2 = Breakpoint(self.__chr_2, self.__pos_2, self.__gene_2, self.__desc_2, self.__logger)  

            if not self.__check_if_bkp_in_panel():
                self.__logger.error("Genes not in panel")
                raise Exception("Genes not in panel")

            if not self.__check_if_bkp_in_coding():
                self.__logger.error("both breakpoints not in coding")
                raise Exception("Both breakpoints not in coding")

            self.__generate_note_variables(__bkp_1,__bkp_2)

        except Exception as e:
            self.__logger.error("Exception in SV annotator")
            self.__logger.error(e)

    def set_is_oncogenic(self,__val): 
        r"""
        This is the public function used by utils to set SV oncokb oncogenic flag
        `__is_oncogenic`

        Parameters
        ----------
        self
        __val : bool

        Returns
        -------
        None
        

        """
        self.__is_oncogenic = __val

    def set_bkpsites(self,__val): 
        r"""
        This is the public function used by utils to set SV breakpoint site

        Parameters
        ----------
        self
        __val : string

        Returns
        -------
        None

        """
        self.__bkpsites = __val

    def __check_if_bkp_in_panel(self):
        r"""
        Check if both breakpoints is in panel

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__partner_1.get_is_panel() or self.__partner_2.get_is_panel()

    def __check_if_bkp_in_coding(self):
        r"""
        Check if any of the breakpoints is in coding

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__partner_1.get_is_panel_and_coding() or self.__partner_2.get_is_panel_and_coding()

    def __check_if_both_breakpoints_in_coding(self):
        r"""
        Check if both the breakpoints are in coding

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__partner_1.get_is_coding() and self.__partner_2.get_is_coding()

    def __check_if_both_breakpoints_in_panel_coding(self):
        r"""
        Check if both the breakpoints are is in panel and in coding

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__partner_1.get_is_panel_and_coding() and self.__partner_2.get_is_panel_and_coding()

    def check_is_intragenic(self): 
        r"""
        Check if sv is intragenic

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return all([self.__partner_1.get_gene() == self.__partner_2.get_gene(), self.__partner_1.get_transcript() == self.__partner_2.get_transcript(), self.__partner_1.get_is_coding(), self.__partner_2.get_is_coding()])
        
    def check_if_sv_type_translocation(self): 
        r"""
        Check if sv type is translocation

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__sv_type == SV_TYPE_TRANSLOCATION

    def check_if_sv_type_deletion(self): 
        r"""
        Check if sv type is deletion

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__sv_type == SV_TYPE_DELETION

    def check_if_sv_type_duplication(self): 
        r"""
        Check if sv type is duplication

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__sv_type == SV_TYPE_DUPLICATION

    def check_if_sv_type_inversion(self): 
        r"""
        Check if sv type is inversion

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__sv_type == SV_TYPE_INVERSION

    def check_if_sv_in_frame(self):
        r"""
        Check if sv is in frame

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return INFRAME_IDENTIFIER in self.__event_info

    def get_is_protein_fusion(self): 
        r"""
        Check if sv is a protien fusion

        Parameters
        ----------
        self

        Returns
        -------
        bool

        Notes
        ------------
        `event_info` should contian `PROTEIN FUSION` or `TRANSCRIPT FUSION`

        """
        return any(x in self.__event_info.upper() for x in [PROTEIN_FUSION_IDENTIFIER,TRANSCRIPT_FUSION_IDENTIFIER]) and  self.__check_if_both_breakpoints_in_coding() 

    def get_is_known_oncokb_fusion(self): 
        r"""
        Check if sv is a known oncokb fusion

        Parameters
        ----------
        self

        Returns
        -------
        bool

        """
        return self.__is_known_oncokb_fusion

    def get_sv_type(self): 
        r"""
        Public method to get the sv type

        Parameters
        ----------
        self

        Returns
        -------
        __sv_type : string

        """
        return self.__sv_type
    
    def get_exons_note(self):
        r"""
        Public method to get exon_note

        Parameters
        ----------
        self

        Returns
        -------
        __exons_note : string

        """
        return self.__exons_note

    def __get_sv_oncokb_type(self): 
        r"""
        This function returns the oncokb sv type `FUSION`, `VII`, `VIII`, `KDD` or `CTD` 

        Parameters
        ----------
        self

        Returns
        -------
        __sv_special_type : string
            oncokb sv type `FUSION`, `VII`, `VIII`, `KDD` or `CTD`

        """
        if self.__fusion_type == NOTE_FUSION_TYPE_FUSION: return SV_SPECICAL_TYPE_FUSION
        elif self.__sv_special_type not in (SV_SPECICAL_TYPE_UNK, SV_SPECICAL_TYPE_FUSION, SV_SPECICAL_TYPE_CDKN2A_16_ISOFORM, SV_SPECICAL_TYPE_CDKN2A_14_ISOFORM, SV_SPECICAL_TYPE_CDKN2A_BOTH_ISOFORM) : 
            return self.__sv_special_type
        return SV_SPECICAL_TYPE_UNK

    def get_annotation_partner_1(self): 
        r"""
        Public method to get annotation partner 1

        Parameters
        ----------
        self

        Returns
        -------
        __partner_1 : Breakpoint

        """
        return self.__partner_1

    def get_annotation_partner_2(self): 
        r"""
        Public method to get annotation partner 2

        Parameters
        ----------
        self

        Returns
        -------
        __partner_2 : Breakpoint

        """
        return self.__partner_2

    def get_annotation(self): 
        r"""
        Public method to get annotation note

        Parameters
        ----------
        self

        Returns
        -------
        annotaion_note : string

        Notes
        ---------
        Combines different parts of the note variable to form annotation

        """
        __fusion_type = "{f}:".format(f=self.__fusion_type) if self.__fusion_type != "" else "" 
        __transcript_coordinates = "{tc}{tc_postfix}{sv_type}".format(tc=self.__transcript_coordinates,tc_postfix=self.__transcript_coordinates_postfix,sv_type=self.__annotation_note_sv_type)
        __note_list = [self.__gene_note, __fusion_type, __transcript_coordinates]
        return self.__join(__note_list)

    def get_comment(self):
        r"""
        Public method to get comment note

        Parameters
        ----------
        self

        Returns
        -------
        __note : string

        Notes
        ---------
        Combines different parts of the note variable to form comment.
        if oncokb sv type in `VIII`, `KDD`, `CDKN2A` or `CTD` return special note from
        constants.py or combine special note with formed note

        """
        __note_list = [COMMENT_PREFIX, self.__comment_gene_note,self.__fusion_type,self.__note_prefix, self.__exons_note, self.__bkpsites, self.__misc_note, self.__sig_note]
        __note = self.__join(__note_list)
        self.__sv_special_type = get_sv_oncokb_type(self)
        if self.__sv_special_type in (SV_SPECICAL_TYPE_VIII,SV_SPECICAL_TYPE_CTD,SV_SPECICAL_TYPE_KDD) : 
            __note = "{prefix} {special_note}".format(prefix=COMMENT_PREFIX_1,special_note=SV_SPECICAL_TYPE_NOTE[self.__sv_special_type])
        elif self.__sv_special_type in (SV_SPECICAL_TYPE_ERG,SV_SPECICAL_TYPE_CDKN2A_BOTH_ISOFORM,SV_SPECICAL_TYPE_CDKN2A_16_ISOFORM,SV_SPECICAL_TYPE_CDKN2A_14_ISOFORM):
            __note = "{note} {special_note}".format(note=__note,special_note=SV_SPECICAL_TYPE_NOTE[self.__sv_special_type])
        return __note if __note.strip() != COMMENT_PREFIX else ""
                           
    def get_position(self):
        r"""
        Public method to get position note

        Parameters
        ----------
        self

        Returns
        -------
        __position : string

        Notes
        ---------
        Combines different parts of the note variable to form position

        """
        if self.get_is_protein_fusion(): 
            __position = "{g1} {exon} {e1} {to} {g2} {exon} {e2}".format(g1=self.__partner_1.get_gene(), g2=self.__partner_2.get_gene(), \
                            e1=self.__partner_1.get_exon(), e2=self.__partner_2.get_exon(),exon=EXONIC_IDENTIFIER.lower(),to=COMMENT_JOIN_STRING_TO)
        else:
            __position = re.sub(r'\.|\bof \b|\binvolves \b|\bwith breakpoints in \b|\bwith a breakpoint in \b', "", self.__exons_note, count=1)
        return __position

    def get_notes(self):
        r"""
        Public method to get annotation, comment, position and oncokb sv type

        Parameters
        ----------
        self

        Returns
        -------
        notes : tuple
            (annotation, comment, position, oncokb_sv_type)

        """
        return (self.get_annotation(),self.get_comment(),self.get_position(),self.__get_sv_oncokb_type())
    

 
if __name__ == '__main__':
    from mskcc.dmp.cvr.utils.dmp_logger import DMPLogger
    logger = DMPLogger()
    logger = logger.getCVRLogger("dmp_sv_load")

    # print SVAnnotator('PCGF3','FAT1','4','4','756508','187518253','Intron of PCGF3(+):1Kb after exon 9','Exon 25 of FAT1(-)','INV','Protein Fusion: mid-exon  {FAT1:PCGF3}',logger).get_notes()
    # print SVAnnotator('TSC1','COL5A1','9','9','135804243','137690052','Exon 3 of TSC1(-)','Intron of COL5A1(+):201bp before exon 37','DEL','-',logger).get_notes()
    # print SVAnnotator('ALK','EML4','2','2','29447038','42493173','Intron of ALK(-):644bp before exon 20','Intron of EML4(+):1Kb after exon 5','INV','Protein Fusion: in frame  {EML4:ALK}',logger).get_notes()
    # print SVAnnotator('EP300','PTPRN2','22','7','41558713','157390026','Intron of EP300(+):13bp before exon 21','Intron of PTPRN2(-):2Kb before exon 16','TRA','Protein Fusion: out of frame  {EP300:PTPRN2}',logger).get_notes()
    # print SVAnnotator('BCR','ABL1','22','9','23632014','133729722','Intron of BCR(+):206bp after exon 13','Intron of ABL1(+):98bp after exon 2','TRA','Protein Fusion: in frame  {BCR:ABL1}',logger).get_notes()
    # print SVAnnotator('KMT2D','KMT2D','12','12','49434833','49441057','Exon 31 of KMT2D(-)','Intron of KMT2D(-):484bp before exon 15','DEL','Deletion within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('ZNF710','IDH2','15','15','90523585','90628575','IGR: 21Kb before ZNF710(+)','Exon 8 of IDH2(-)','DEL','-',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39798230','42869522','Intron of ERG(-):3Kb before exon 3','Intron of TMPRSS2(-):523bp after exon 2','DEL','Protein Fusion: out of frame  {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('NOTCH3','PSMC5','19','17','15272285','61904863','Exon 33 of NOTCH3(-)','Exon 1 of PSMC5(+)','TRA','Protein Fusion: mid-exon  {NOTCH3:PSMC5}',logger).get_notes()
    # print SVAnnotator('CTNNA3','FGFR2','10','10','68009575','123240911','Intron of CTNNA3(-):31Kb after exon 13','Intron of FGFR2(-):1Kb before exon 18','DEL','Protein Fusion: in frame  {FGFR2:CTNNA3}',logger).get_notes()
    # print SVAnnotator('ZC3H7B','UPF1','22','19','41744573','18976550','Intron of ZC3H7B(+):402bp after exon 15','Exon 22 of UPF1(+)','TRA','Protein Fusion: mid-exon  {UPF1:ZC3H7B}',logger).get_notes()
    # print SVAnnotator('GATA1','WNK3','X','X','48650328','54291385','Exon 3 of GATA1(+)','Intron of WNK3(-):6Kb before exon 11','DUP','-',logger).get_notes()
    # print SVAnnotator('KMT2D','FAM149B1','12','10','49446311','74969560','Intron of KMT2D(-):35bp after exon 9','Intron of FAM149B1(+):448bp before exon 7','TRA','Protein Fusion: out of frame  {KMT2D:FAM149B1}',logger).get_notes()
    # print SVAnnotator('CIC','PAFAH1B3','19','19','42793198','42806688','Exon 7 of CIC(+)','5-UTR of PAFAH1B3(-):5Kb before coding start','DUP','-',logger).get_notes()
    # print SVAnnotator('RASA1','RASA1','5','5','86574220','86659309','Intron of RASA1(+):9Kb after exon 1','Exon 11 of RASA1(+)','DUP','Duplication within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('ARID1A','GPATCH3','1','1','27102127','27222129','Exon 19 of ARID1A(+)','Intron of GPATCH3(-):1Kb before exon 3','INV','Protein Fusion: mid-exon  {ARID1A:GPATCH3}',logger).get_notes()
    # print SVAnnotator('ZNF316','ETV1','7','7','6669540','14025715','IGR: 7Kb before ZNF316(+)','Intron of ETV1(-):547bp after exon 4','DEL','-',logger).get_notes()
    # print SVAnnotator('ZFHX3','ZFHX3','16','16','72798880','72829520','IGR: 18Kb before ZFHX3(-)','Exon 9 of ZFHX3(-)','DEL','-',logger).get_notes()
    # print SVAnnotator('PAK1','MIR7976','11','11','77046959','97205904','Intron of PAK1(-):171bp after exon 13','IGR: 285Kb before MIR7976(-)','DEL','-',logger).get_notes()
    # print SVAnnotator('ATXN7L3B','BIRC3','12','11','74859022','102195780','IGR: 73Kb before ATXN7L3B(+)','Exon 3 of BIRC3(+)','TRA','-',logger).get_notes()
    # print SVAnnotator('PTCHD1','NAB2','X','12','23509870','57486173','IGR: 157Kb before PTCHD1(+)','Intron of NAB2(+):57bp before exon 3','TRA','-',logger).get_notes()
    # print SVAnnotator('EIF1AY','TBX3','Y','12','22716517','115109879','IGR: 21Kb before EIF1AY(+)','Exon 8 of TBX3(-)','TRA','-',logger).get_notes()
    # print SVAnnotator('TBX3','TBX3','12','12','115111461','115119062','Intron of TBX3(-):508bp after exon 7','Intron of TBX3(-):111bp before exon 2','INV','-',logger).get_notes()
    # print SVAnnotator('TFE3','ASPSCR1','X','17','48893253','79967177','Intron of TFE3(-):1Kb before exon 6','Intron of ASPSCR1(+):110bp after exon 8','TRA','Protein Fusion: out of frame  {TFE3:ASPSCR1}',logger).get_notes()
    # print SVAnnotator('STK19','C4A','6','6','31948448','31963776','Exon 7 of STK19(+)','Exon 26 of C4A(+)','DEL','Protein Fusion: mid-exon  {STK19:C4A}',logger).get_notes()
    # print SVAnnotator('PIK3C2G','FBXO16','12','8','18552829','28324745','Intron of PIK3C2G(+):48bp after exon 15','Intron of FBXO16(-):3Kb before exon 3','TRA','Protein Fusion: out of frame  {PIK3C2G:FBXO16}',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39871604','42871198','Promoter of ERG(-):120Kb from tx start','Intron of TMPRSS2(-):1Kb before exon 2','DEL','Transcript Fusion {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('ELF3','ELF3','1','1','201979985','201982188','5-UTR of ELF3(+):279bp before coding start','Intron of ELF3(+):24bp after exon 6','INV','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('ROS1','ROS1','6','6','117725507','117731222','Exon 5 of ROS1(-)','Intron of ROS1(-):417bp before exon 4','INV','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('TSC2','TSC2','16','16','2134097','2139229','Intron of TSC2(+):131bp before exon 34','Promoter of TSC2(+):41Kb from tx start','DEL','Deletion within transcript',logger).get_notes()
    # print SVAnnotator('KMT2B','KMT2E','19','7','36224589','104695690','Intron of KMT2B(+):2bp after exon 29','Intron of KMT2E(+):7Kb before exon 3','TRA','Protein Fusion: in frame  {KMT2B:KMT2E}',logger).get_notes()
    # print SVAnnotator('CA10','RECQL4','17','8','50212204','145742559','Intron of CA10(-):23Kb after exon 1','Exon 4 of RECQL4(-)','TRA','Protein Fusion: mid-exon  {CA10:RECQL4}',logger).get_notes()
    # print SVAnnotator('EPHA5','EPHA5','4','4','66241266','66242724','Intron of EPHA5(-):1Kb after exon 9','Exon 9 of EPHA5(-)','DEL','Deletion within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('KMT2C','KMT2C','7','7','151827709','151859482','IGR: 4Kb before KMT2C(-)','Exon 43 of KMT2C(-)','DEL','-',logger).get_notes()
    # print SVAnnotator('EGFR','GBAS','7','7','55268444','56048007','Intron of EGFR(+):338bp after exon 24','Intron of GBAS(+):1Kb before exon 4','DUP','Protein Fusion: out of frame  {GBAS:EGFR}',logger).get_notes()
    # print SVAnnotator('TP53BP1','RNU6-28P','15','15','43785386','43810133','Promoter of TP53BP1(-):86Kb from tx start','5-UTR of RNU6-28P(+):85Kb before coding start','DEL','-',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39941388','42872642','Intron of ERG(-):6Kb after exon 3','Intron of TMPRSS2(-):3Kb before exon 2','DEL','Protein Fusion: out of frame  {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('MET','ST7','7','7','116399527','116756324','Exon 10 of MET(+)','Intron of ST7(+):3Kb before exon 3','DUP','Protein Fusion: mid-exon  {ST7:MET}',logger).get_notes()
    # print SVAnnotator('TJP1','SPRED1','15','15','30078906','38456891','Intron of TJP1(-):13Kb before exon 4','IGR: 88Kb before SPRED1(+)','DUP','-',logger).get_notes()
    # print SVAnnotator('POLB','AGO2','8','8','42201644','141554151','Intron of POLB(+):826bp before exon 3','Intron of AGO2(-):160bp after exon 14','INV','Protein Fusion: out of frame  {AGO2:POLB}',logger).get_notes()
    # print SVAnnotator('NF1','TOM1L1','17','17','29557072','52929124','Intron of NF1(+):80bp after exon 22','IGR: 49Kb before TOM1L1(+)','INV','-',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39861145','42876260','Intron of ERG(-):9Kb after exon 1','Intron of TMPRSS2(-):4Kb after exon 1','DEL','Protein Fusion: out of frame  {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('PDCD1','RTP5','2','2','242793230','242810759','Exon 5 of PDCD1(-)','Promoter of RTP5(+):1Kb from tx start','INV','Protein Fusion: mid-exon  {RTP5:PDCD1}',logger).get_notes()
    # print SVAnnotator('CDKN2A','CDKN2B-AS1','9','9','21974509','22066617','Exon 1 of CDKN2A(-)','5-UTR of CDKN2B-AS1(+):54Kb before coding start','DEL','-',logger).get_notes()
    # print SVAnnotator('CTNNB1','CTNNB1','3','3','41265013','41266391','5-UTR of CTNNB1(+):546bp before coding start','Intron of CTNNB1(+):53bp before exon 4','DEL','Deletion within transcript',logger).get_notes()
    # print SVAnnotator('PIK3R3','PIK3R3','1','1','46531693','46563923','Intron of PIK3R3(-):32bp after exon 5','Intron of PIK3R3(-):18Kb before exon 2','DUP','Duplication of 4 exons : out of frame',logger).get_notes()
    # print SVAnnotator('CDKN2A','CDKN2A','9','9','21968241','21970901','Exon 3 of CDKN2A(-)','Exon 2 of CDKN2A(-)','DEL','Deletion within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('STAC2','STAT5B','17','17','37367870','40359544','3-UTR of STAC2(-):674bp after coding stop','Intron of STAT5B(-):31bp after exon 16','INV','-',logger).get_notes()
    # print SVAnnotator('ANKRD7','TAP2','7','6','118149779','32798010','IGR: 285Kb before ANKRD7(+)','Intron of TAP2(-):33bp after exon 9','TRA','-',logger).get_notes()
    # print SVAnnotator('SLX4','CREBBP','16','16','3644966','3823744','Intron of SLX4(-):366bp before exon 10','Intron of CREBBP(-):7bp after exon 13','DEL','Protein Fusion: in frame  {CREBBP:SLX4}',logger).get_notes()
    # print SVAnnotator('DNMT1','ZNF701','19','19','10267093','53065226','Exon 17 of DNMT1(-)','IGR: 8Kb before ZNF701(+)','DEL','-',logger).get_notes()
    # print SVAnnotator('KDM6A','KDM6A','X','X','44761938','44950141','Intron of KDM6A(+):29Kb after exon 2','Intron of KDM6A(+):32bp after exon 26','DUP','Duplication of 24 exons : out of frame',logger).get_notes()
    # print SVAnnotator('NPAS1','BBC3','19','19','47535039','47727316','Intron of NPAS1(+):496bp before exon 3','Intron of BBC3(-):2Kb before exon 4','INV','Protein Fusion: in frame  {NPAS1:BBC3}',logger).get_notes()
    # print SVAnnotator('NTRK3','ETV6','15','12','88522437','12037217','Intron of NTRK3(-):38Kb before exon 15','Intron of ETV6(+):161bp before exon 6','TRA','Protein Fusion: in frame  {ETV6:NTRK3}',logger).get_notes()
    # print SVAnnotator('KMT2C','KMT2C','7','7','151875002','151877411','Exon 38 of KMT2C(-)','Intron of KMT2C(-):200bp before exon 37','DEL','Deletion within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('BCR','ABL1','22','9','23634581','133586273','Intron of BCR(+):146bp before exon 15','Promoter of ABL1(+):3Kb from tx start','TRA','Transcript Fusion {ABL1:BCR}',logger).get_notes()
    # print SVAnnotator('S100A10','NTRK1','1','1','151930587','156835291','IGR: 25Kb before S100A10(-)','Intron of NTRK1(+):700bp after exon 3','DEL','-',logger).get_notes()
    # print SVAnnotator('PIK3R3','LRRC41','1','1','46509488','46756288','Exon 10 of PIK3R3(-)','Intron of LRRC41(-):4Kb before exon 4','DUP','Protein Fusion: mid-exon  {PIK3R3:LRRC41}',logger).get_notes()
    # print SVAnnotator('MYC','SNRK','8','3','128750517','43357294','Exon 2 of MYC(+)','Intron of SNRK(+):12Kb after exon 3','TRA','Protein Fusion: mid-exon  {SNRK:MYC}',logger).get_notes()
    # print SVAnnotator('IGF1R','MAP3K13','15','3','99500351','185040460','Exon 21 of IGF1R(+)','Intron of MAP3K13(+):37Kb after exon 2','TRA','Protein Fusion: mid-exon  {MAP3K13:IGF1R}',logger).get_notes()
    # print SVAnnotator('FAM175A','CENPE','4','4','84390340','104032264','Intron of FAM175A(-):36bp before exon 6','Intron of CENPE(-):96bp before exon 45','DEL','Protein Fusion: out of frame  {CENPE:FAM175A}',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39861753','42870567','Intron of ERG(-):9Kb after exon 1','Intron of TMPRSS2(-):451bp before exon 2','DEL','Protein Fusion: out of frame  {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('RTEL1','RTEL1','20','20','62310376','62320863','Intron of RTEL1(+):677bp after exon 12','Exon 23 of RTEL1(+)','DUP','Duplication within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('EPS15','NTRK1','1','1','51858813','156844212','Intron of EPS15(-):1Kb after exon 21','Intron of NTRK1(+):20bp after exon 9','INV','Protein Fusion: in frame  {EPS15:NTRK1}',logger).get_notes()
    # print SVAnnotator('MIR5702','IRS1','2','2','227329627','227662807','IGR: 194Kb before MIR5702(-)','Exon 1 of IRS1(-)','DUP','-',logger).get_notes()
    # print SVAnnotator('HS6ST3','HIST1H1C','13','6','97043138','26056197','Intron of HS6ST3(+):299Kb after exon 1','Exon 1 of HIST1H1C(-)','TRA','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('RTEL1','ZBTB46','20','20','62326916','62425706','Intron of RTEL1(+):83bp after exon 34','5-UTR of ZBTB46(-):47Kb before coding start','DEL','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('STAG2','STAG2','X','X','123108022','123156541','5-UTR of STAG2(+):48Kb before coding start','Intron of STAG2(+):20bp after exon 3','INV','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('SOCS1','SOCS1','16','16','11348944','11350617','Exon 2 of SOCS1(-)','Promoter of SOCS1(-):2Kb from tx start','INV','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('POLE','SLC9A2','12','2','133243918','103227722','Intron of POLE(-):170bp after exon 20','IGR: 8Kb before SLC9A2(+)','TRA','-',logger).get_notes()
    # print SVAnnotator('CDH1','CDH1','16','16','68770959','68772126','Promoter of CDH1(+):235bp from tx start','Intron of CDH1(+):73bp before exon 2','DEL','Deletion within transcript',logger).get_notes()
    # print SVAnnotator('MYD88','OXSR1','3','3','38182163','38214127','Intron of MYD88(+):84bp before exon 4','Intron of OXSR1(+):7Kb after exon 1','DUP','Protein Fusion: out of frame  {OXSR1:MYD88}',logger).get_notes()
    # print SVAnnotator('ALK','EML4','2','2','29446522','42496505','Intron of ALK(-):128bp before exon 20','Intron of EML4(+):5Kb after exon 5','INV','Protein Fusion: in frame  {EML4:ALK}',logger).get_notes()
    # print SVAnnotator('TAP1','HLA-DMB','6','6','32815732','32904545','Exon 8 of TAP1(-)','Intron of HLA-DMB(-):403bp after exon 3','DEL','Protein Fusion: mid-exon  {HLA-DMB:TAP1}',logger).get_notes()
    # print SVAnnotator('FER1L4','PTPRT','20','20','34150608','40710478','3-UTR of FER1L4(-):45Kb after coding stop','Intron of PTPRT(-):43bp after exon 31','DEL','-',logger).get_notes()
    # print SVAnnotator('TNFRSF14','PLXNA2','1','1','2489056','208364734','Intron of TNFRSF14(+):108bp before exon 2','Intron of PLXNA2(-):19Kb after exon 3','INV','Protein Fusion: in frame  {TNFRSF14:PLXNA2}',logger).get_notes()
    # print SVAnnotator('ETV6','BCL2L14','12','12','12006418','12150359','Exon 4 of ETV6(+)','IGR: 74Kb before BCL2L14(+)','DUP','-',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55134073','55222656','Intron of EGFR(+):47Kb after exon 1','Intron of EGFR(+):811bp after exon 7','DEL','Deletion of 6 exons : in frame',logger).get_notes()
    # print SVAnnotator('PRKACA','DNAJB1','19','19','14219638','14627553','Intron of PRKACA(-):1Kb before exon 2','Exon 2 of DNAJB1(-)','DEL','Protein Fusion: mid-exon  {DNAJB1:PRKACA}',logger).get_notes()
    # print SVAnnotator('CASP8','CASP8','2','2','202108480','202141637','5-UTR of CASP8(+):23Kb before coding start','Exon 7 of CASP8(+)','DEL','Deletion within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('APC','SRP19','5','5','112176002','112200782','Exon 16 of APC(+)','Intron of SRP19(+):353bp after exon 4','DEL','Protein Fusion: mid-exon  {APC:SRP19}',logger).get_notes()
    # print SVAnnotator('ALK','EML4','2','2','29446646','42527742','Intron of ALK(-):252bp before exon 20','Intron of EML4(+):638bp before exon 13','INV','Protein Fusion: in frame  {EML4:ALK}',logger).get_notes()
    # print SVAnnotator('CDH1','CDH1','16','16','68848903','68855972','Intron of CDH1(+):514bp before exon 10','Exon 12 of CDH1(+)','DUP','Duplication within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('RET','TIMM23B','10','10','43611448','51584256','Intron of RET(+):583bp before exon 12','Intron of TIMM23B(+):149Kb before exon 7','DUP','Protein Fusion: out of frame  {TIMM23B:RET}',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39911122','42868304','Intron of ERG(-):36Kb after exon 3','Intron of TMPRSS2(-):2Kb after exon 2','DEL','Protein Fusion: in frame  {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('NEGR1','NEGR1','1','1','72242069','72554634','Intron of NEGR1(-):89bp before exon 3','Intron of NEGR1(-):154Kb before exon 2','DEL','Deletion of 1 exon : out of frame',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39863298','42869384','Intron of ERG(-):7Kb after exon 1','Intron of TMPRSS2(-):661bp after exon 2','DEL','Protein Fusion: in frame  {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('RPTOR','RPTOR','17','17','78919515','78926053','Exon 26 of RPTOR(+)','Intron of RPTOR(+):3Kb after exon 28','DEL','Deletion within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('RICTOR','TTC33','5','5','38946536','40729443','Intron of RICTOR(-):33bp after exon 33','Intron of TTC33(-):865bp before exon 4','DUP','Protein Fusion: out of frame  {RICTOR:TTC33}',logger).get_notes()
    # print SVAnnotator('SMARCA2','CD274','9','9','2174635','5457097','Intron of SMARCA2(+):4Kb after exon 29','Exon 3 of CD274(+)','INV','-',logger).get_notes()
    # print SVAnnotator('TERT','FAM173B','5','5','1295168','10249820','Promoter of TERT(-):42Kb from tx start','Intron of FAM173B(-):149bp after exon 1','DEL','Transcript Fusion {FAM173B:TERT}',logger).get_notes()
    # print SVAnnotator('KMT2C','KMT2C','7','7','151854877','151952175','Intron of KMT2C(-):1Kb after exon 44','Intron of KMT2C(-):2Kb before exon 10','INV','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('TMPRSS2','TMPRSS2','21','21','42857943','42870416','Intron of TMPRSS2(-):2Kb after exon 5','Intron of TMPRSS2(-):300bp before exon 2','INV','Antisense Fusion',logger).get_notes()
    # print SVAnnotator('ERG','TMPRSS2','21','21','39936477','42870000','Intron of ERG(-):11Kb after exon 3','Intron of TMPRSS2(-):45bp after exon 2','DEL','Protein Fusion: in frame  {TMPRSS2:ERG}',logger).get_notes()
    # print SVAnnotator('BNC2','TEK','9','9','16340059','27169375','IGR: 69Kb before BNC2(-)','Intron of TEK(+):99bp before exon 4','INV','-',logger).get_notes()
    # print SVAnnotator('SLC35E3','MDM2','12','12','69140962','69233041','Intron of SLC35E3(+):403bp after exon 1','Intron of MDM2(+):12bp before exon 11','DEL','Protein Fusion: in frame  {SLC35E3:MDM2}',logger).get_notes()
    # print SVAnnotator('C12orf79','BTG1','12','12','92278193','92538071','IGR: 101Kb before C12orf79(-)','Exon 2 of BTG1(-)','DEL','-',logger).get_notes()
    # print SVAnnotator('BCR','ABL1','22','9','23632014','133729722','Intron of BCR(+):206bp after exon 13','Intron of ABL1(+):98bp after exon 2','TRA','Protein Fusion: in frame  {BCR:ABL1}',logger).get_notes()
    # print SVAnnotator('CDK12','FGFR4','17','5','37687238','176518227','Exon 14 of CDK12(+)','Intron of FGFR4(+):122bp after exon 5','TRA','Protein Fusion: mid-exon  {FGFR4:CDK12}',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55241243','55269327','Intron of EGFR(+):370bp before exon 18','Intron of EGFR(+):100bp before exon 26','DUP','Duplication of 8 exons',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55268516','55269654','Intron of EGFR(+):364bp before exon 25','Intron of EGFR(+):179bp after exon 26','DEL','Deletion of 2 exons : in frame',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55268265','55269654','Intron of EGFR(+):159bp after exon 24','Intron of EGFR(+):179bp after exon 26','DEL','Deletion of 2 exons : in frame',logger).get_notes()
    # print SVAnnotator('KIAA1549','BRAF','7','7','138548248','140491093','Intron of KIAA1549(-):2Kb before exon 16','Intron of BRAF(-):3Kb after exon 8','DUP','Protein Fusion: in frame  {KIAA1549:BRAF}',logger).get_notes()
    # print SVAnnotator('ALK','EML4','2','2','29447383','42503133','Intron of ALK(-):943bp after exon 19','Intron of EML4(+):5Kb before exon 6','INV','Protein Fusion: in frame  {EML4:ALK}',logger).get_notes()
    # print SVAnnotator('ALK','EML4','2','2','29446514','42553077','Intron of ALK(-):120bp before exon 20','Intron of EML4(+):216bp before exon 20','INV','Protein Fusion: in frame  {EML4:ALK}',logger).get_notes()
    # print SVAnnotator('NTRK1','ATP1A2','1','1','156844714','160106162','Exon 11 of NTRK1(+)','Intron of ATP1A2(+):2bp after exon 18','DUP','Protein Fusion: mid-exon  {ATP1A2:NTRK1}',logger).get_notes()
    # print SVAnnotator('NTRK3','ETV6','15','12','88651251','12020963','Intron of NTRK3(-):18Kb after exon 13','Intron of ETV6(+):1Kb before exon 5','TRA','Protein Fusion: in frame  {ETV6:NTRK3}',logger).get_notes()
    # print SVAnnotator('PRKACA','DNAJB1','19','19','14227329','14627970','Intron of PRKACA(-):984bp after exon 1','Intron of DNAJB1(-):112bp before exon 2','DEL','Protein Fusion: in frame  {DNAJB1:PRKACA}',logger).get_notes()
    # print SVAnnotator('WDR11','FGFR2','10','10','122913765','123242665','IGR: 303Kb before WDR11(+)','Intron of FGFR2(-):546bp after exon 17','INV','-',logger).get_notes()
    # print SVAnnotator('NTRK2','CENPP','9','9','87563367','95141254','Intron of NTRK2(+):9bp before exon 18','Intron of CENPP(+):790bp before exon 4','DUP','Protein Fusion: out of frame  {CENPP:NTRK2}',logger).get_notes()
    # print SVAnnotator('MKRN1','BRAF','7','7','140157974','140488436','Intron of MKRN1(-):832bp after exon 4','Intron of BRAF(-):1Kb before exon 9','DUP','Protein Fusion: in frame  {MKRN1:BRAF}',logger).get_notes()
    # print SVAnnotator('EEA1','EGFR','12','7','93221108','55241575','Intron of EEA1(-):579bp after exon 12','Intron of EGFR(+):38bp before exon 18','TRA','Protein Fusion: in frame  {EEA1:EGFR}',logger).get_notes()
    # print SVAnnotator('ROS1','DCBLD1','6','6','117643817','117885282','Intron of ROS1(-):1Kb before exon 35','Intron of DCBLD1(+):6Kb before exon 15','DEL','-',logger).get_notes()
    # print SVAnnotator('NTRK3','ETV6','15','12','88505787','12029143','Intron of NTRK3(-):22Kb before exon 15','Intron of ETV6(+):6Kb after exon 5','TRA','Protein Fusion: in frame  {ETV6:NTRK3}',logger).get_notes()
    # print SVAnnotator('NTRK1','PLEKHA6','1','1','156844239','204214638','Intron of NTRK1(+):47bp after exon 9','Intron of PLEKHA6(-):104bp after exon 14','INV','Protein Fusion: in frame  {PLEKHA6:NTRK1}',logger).get_notes()
    # print SVAnnotator('MAP3K13','MASP1','3','3','185161302','187005422','Exon 4 of MAP3K13(+)','Intron of MASP1(-):2Kb before exon 2','INV','Protein Fusion: mid-exon  {MAP3K13:MASP1}',logger).get_notes()
    # print SVAnnotator('NUMA1','EGFR','11','7','71728415','55241128','Intron of NUMA1(-):317bp after exon 13','Intron of EGFR(+):311bp after exon 17','TRA','Protein Fusion: in frame  {NUMA1:EGFR}',logger).get_notes()
    # print SVAnnotator('ALK','EML4','2','2','29446930','42526326','Intron of ALK(-):536bp before exon 20','Intron of EML4(+):2Kb before exon 13','INV','Protein Fusion: in frame  {EML4:ALK}',logger).get_notes()
    # print SVAnnotator('ALK','EML4','2','2','29447857','42522763','Intron of ALK(-):469bp after exon 19','Intron of EML4(+):107bp after exon 12','INV','Protein Fusion: in frame  {EML4:ALK}',logger).get_notes()
    # print SVAnnotator('ROS1','DCBLD1','6','6','117644742','117886239','Intron of ROS1(-):752bp after exon 34','Intron of DCBLD1(+):5Kb before exon 15','DEL','-',logger).get_notes()
    # print SVAnnotator('TACC3','FGFR3','4','4','1725566','1808902','Exon 3 of TACC3(+)','Exon 18 of FGFR3(+)','DUP','Protein Fusion: mid-exon  {FGFR3:TACC3}',logger).get_notes()
    # print SVAnnotator('FGFR2','TACC2','10','10','123239743','123975434','Intron of FGFR2(-):208bp before exon 18','Intron of TACC2(+):468bp after exon 9','INV','Protein Fusion: in frame  {FGFR2:TACC2}',logger).get_notes()
    # print  SVAnnotator('PPP1R9A','MET','7','7','94865834','116422185','Intron of PPP1R9A(+):10Kb after exon 6','Intron of MET(+):34bp after exon 18','DUP','Protein Fusion: out of frame  {MET:PPP1R9A}',logger).get_notes()
    # print SVAnnotator('RET','CCDC6','10','10','43611418','61639380','Intron of RET(+):613bp before exon 12','Intron of CCDC6(-):26Kb after exon 1','INV','Protein Fusion: in frame  {CCDC6:RET}',logger).get_notes()
    # print SVAnnotator('KIF5B','RET','10','10','32316218','43610413','Intron of KIF5B(-):1Kb after exon 15','Intron of RET(+):229bp after exon 11','INV','Protein Fusion: in frame  {KIF5B:RET}',logger).get_notes()
    # print  SVAnnotator('ERBB2','GSDMA','17','17','37865415','38125047','Intron of ERBB2(+):155bp before exon 4','Intron of GSDMA(+):2Kb before exon 4','DUP','Protein Fusion: out of frame  {GSDMA:ERBB2}',logger).get_notes()
    # print SVAnnotator('SND1','BRAF','7','7','127425525','140493552','Intron of SND1(+):22Kb before exon 11','Intron of BRAF(-):555bp after exon 8','INV','Protein Fusion: in frame  {SND1:BRAF}',logger).get_notes()
    # print SVAnnotator('TACC3','FGFR3','4','4','1727442','1804890','Intron of TACC3(+):2Kb after exon 3','Intron of FGFR3(+):528bp before exon 8','DUP','Protein Fusion: out of frame  {FGFR3:TACC3}',logger).get_notes()
    # print SVAnnotator('KIF5B','RET','10','10','32313695','43611965','Intron of KIF5B(-):2Kb before exon 16','Intron of RET(+):66bp before exon 12','INV','Protein Fusion: in frame  {KIF5B:RET}',logger).get_notes()
    # print SVAnnotator('SPECC1L','RET','22','10','24739056','43610353','Intron of SPECC1L(+):4Kb before exon 11','Intron of RET(+):169bp after exon 11','TRA','Protein Fusion: in frame  {SPECC1L:RET}',logger).get_notes()
    # print SVAnnotator('BRAF','KLHL12','7','1','140493984','202880970','Intron of BRAF(-):123bp after exon 8','Intron of KLHL12(-):639bp before exon 5','TRA','Protein Fusion: in frame  {KLHL12:BRAF}',logger).get_notes()
    # print SVAnnotator('KIAA1549','BRAF','7','7','138546697','140489537','Intron of KIAA1549(-):495bp before exon 16','Intron of BRAF(-):2Kb before exon 9','DUP','Protein Fusion: in frame  {KIAA1549:BRAF}',logger).get_notes()
    # print SVAnnotator('RARA','MAP3K3','17','17','38498999','61752118','Intron of RARA(+):6Kb before exon 3','Intron of MAP3K3(+):7Kb before exon 7','DEL','Protein Fusion: in frame  {RARA:MAP3K3}',logger).get_notes()
    # print SVAnnotator('MAP2K2','MAP2K2','19','19','4120570','4123788','Intron of MAP2K2(-):3Kb before exon 2','Exon 1 of MAP2K2(-)','DEL','Deletion within transcript : mid-exon',logger).get_notes()
    # print SVAnnotator('ROS1','CD74','6','5','117646032','149782990','Intron of ROS1(-):454bp before exon 34','Intron of CD74(-):115bp before exon 7','TRA','Protein Fusion: in frame  {CD74:ROS1}',logger).get_notes()
    # print SVAnnotator('TACC3','FGFR3','4','4','1737170','1808872','Intron of TACC3(+):118bp after exon 7','Exon 18 of FGFR3(+)','DUP','Protein Fusion: mid-exon  {FGFR3:TACC3}',logger).get_notes()
    # print SVAnnotator('ROS1','TFG','6','3','117643761','100448554','Intron of ROS1(-):1Kb before exon 35','Intron of TFG(+):852bp after exon 4','TRA','Protein Fusion: in frame  {TFG:ROS1}',logger).get_notes()
    # print SVAnnotator('AURKA','FAM114A2','20','5','54945351','153381607','Exon 9 of AURKA(-)','Intron of FAM114A2(-):203bp after exon 11','TRA','Protein Fusion: mid-exon  {AURKA:FAM114A2}',logger).get_notes()
    # print SVAnnotator('MAP2K4','SLC47A1','17','17','12013590','19443565','Intron of MAP2K4(+):101bp before exon 6','Intron of SLC47A1(+):2Kb before exon 2','DUP','Protein Fusion: in frame  {SLC47A1:MAP2K4}',logger).get_notes()
    # print SVAnnotator('FADD','EGFR','11','7','70052330','55241077','Exon 2 of FADD(+)','Intron of EGFR(+):260bp after exon 17','TRA','Protein Fusion: mid-exon  {EGFR:FADD}',logger).get_notes()
    # print SVAnnotator('KIF5B','RET','10','10','32315869','43612025','Intron of KIF5B(-):1Kb after exon 15','Intron of RET(+):6bp before exon 12','INV','Protein Fusion: in frame  {KIF5B:RET}',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55209999','55221934','Exon 2 of EGFR(+)','Intron of EGFR(+):89bp after exon 7','DUP','Duplication within transcript : mid-exon',logger).get_notes()
    # print  SVAnnotator('BRAF','PRIM2','7','6','140493232','57410815','Intron of BRAF(-):875bp after exon 8','Intron of PRIM2(+):12Kb after exon 10','TRA','Protein Fusion: in frame  {PRIM2:BRAF}',logger).get_notes()
    # print SVAnnotator('NTRK3','ETV6','15','12','88666269','12015303','Intron of NTRK3(-):3Kb after exon 13','Intron of ETV6(+):7Kb before exon 5','TRA','Protein Fusion: in frame  {ETV6:NTRK3}',logger).get_notes()
    # print SVAnnotator('EXOC4','BRAF','7','7','133408363','140490241','Intron of EXOC4(+):93Kb after exon 10','Intron of BRAF(-):3Kb before exon 9','INV','Protein Fusion: out of frame  {EXOC4:BRAF}',logger).get_notes()
    # print SVAnnotator('BRAF','AGK','7','7','140494309','141266226','Intron of BRAF(-):42bp before exon 8','Intron of AGK(+):11Kb after exon 2','INV','Protein Fusion: in frame  {AGK:BRAF}',logger).get_notes()
    # print SVAnnotator('LACE1','ROS1','6','6','108749238','117660863','Intron of LACE1(+):19Kb before exon 8','Intron of ROS1(-):1Kb after exon 30','INV','Protein Fusion: in frame  {LACE1:ROS1}',logger).get_notes()
    # print SVAnnotator('SYT16','AKT1','14','14','62562584','105246431','Intron of SYT16(+):5Kb before exon 6','Exon 3 of AKT1(-)','INV','Protein Fusion: mid-exon  {SYT16:AKT1}',logger).get_notes()
    # print SVAnnotator('TACC3','FGFR3','4','4','1739309','1808878','Intron of TACC3(+):15bp before exon 10','Exon 18 of FGFR3(+)','DUP','Protein Fusion: mid-exon  {FGFR3:TACC3}',logger).get_notes()
    # print SVAnnotator('MAPK3','MAPK3','16','16','30128466','30130544','Intron of MAPK3(-):8bp after exon 6','Intron of MAPK3(-):685bp before exon 3','DUP','Duplication of 4 exons : out of frame',logger).get_notes()
    # print SVAnnotator('CDK12','ERBB2','17','17','37673710','37881361','Exon 10 of CDK12(+)','Exon 21 of ERBB2(+)','DUP','Protein Fusion: mid-exon  {ERBB2:CDK12}',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55207931','55221978','Intron of EGFR(+):2Kb before exon 2','Intron of EGFR(+):133bp after exon 7','DEL','Deletion of 6 exons : in frame',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55197558','55221881','Intron of EGFR(+):12Kb before exon 2','Intron of EGFR(+):36bp after exon 7','DEL','Deletion of 6 exons : in frame',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55209721','55223307','Intron of EGFR(+):257bp before exon 2','Intron of EGFR(+):215bp before exon 8','DEL','Deletion of 6 exons : in frame',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55118037','55223239','Intron of EGFR(+):31Kb after exon 1','Intron of EGFR(+):283bp before exon 8','DEL','Deletion of 6 exons : in frame',logger).get_notes()
    # print SVAnnotator('EGFR','EGFR','7','7','55269624','55272030','Intron of EGFR(+):149bp after exon 26','Intron of EGFR(+):918bp before exon 28','DEL','Deletion of 1 exon : out of frame',logger).get_notes()
    # print SVAnnotator('EGFR','LOC650226','7','7','55268790','56436684','Intron of EGFR(+):90bp before exon 25','IGR: 55Kb before LOC650226(-)','INV','-',logger).get_notes()