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

from constants.constants import *
from utils import *

class Breakpoint(object):
    """
    Class to represent a breakpoint and other related attributes and features.
    """
    __is_entire_kinase = 0
    __is_tumour_suppressor = False
    __first_exon = 1
    __last_exon = ''
    __start_pos = ''
    __stop_pos = ''
    __exon = None
    __intron = None
    __site = None
    __variant_site_1 = None
    __variant_site_2 = None
    __transcript = []

    def __init__(self, chrom, pos, gene, desc, logger):
        r"""
        Constrcutor for Breakpoint 

        Parameters
        ----------
        self
        chrom : string
        pos : int
        gene : string
        desc : string
        logger : logger

        Returns
        -------
        None

        Raises
        -------
        Exception
            Critical exception while getting create attributes for Breakpoint class

        """
        try:
            self.__chrom = str(chrom)
            self.__pos = int(pos)
            self.__gene = str(gene)
            self.__desc = str(desc)
            self.__logger = logger
            self.__expand()
        except TypeError as te:
            logger.error("Could not create a new instance of bkp class due to inappropriate values for parameters.")
            logger.error(te)
            raise ("Could not create a new instance of bkp class")
        except Exception as e:
            logger.error("Unexpected error: {}".format(e))
            raise ("Could not create a new instance of bkp class")

    def __set_strand(self):
        r"""
        Sets strand info based on desc

        Parameters
        ----------
        self

        Returns
        -------
        None

        Raises
        -------
        Exception
            Critical exception: Cannot find strand info in desc

        """
        if POSITIVE_STRAND_IDENTIFIER in self.__desc:
            self.__strand = POSITIVE_STRAND_IDENTIFIER
        elif NEGATIVE_STRAND_IDENTIFIER in self.__desc:
            self.__strand = NEGATIVE_STRAND_IDENTIFIER
        else:
            raise Exception("Cannot get strand info for breakpoint %s:%s."
                            % (self.__chrom, self.__pos))

    def __set_is_coding(self):
        r"""
        Sets __is_coding flag based on cdna and transcript 

        Parameters
        ----------
        self

        Returns
        -------
        None

        """
        self.__is_coding = True if self.__cdna and self.__cdna.startswith("c.") and self.__transcript != "" else False

    def __set_is_panel(self):
        r"""
        Sets __is_panel flag based on if gene is in IMAPCTV6_PANEL in constants file

        Parameters
        ----------
        self

        Returns
        -------
        None

        """
        self.__is_panel = True if self.__gene in IMAPCTV6_PANEL else False

    def __set_is_hotspot(self):
        r"""
        Sets __is_hotspot flag based on if gene is in IMPACTV6_HOTSPOTS in constants file

        Parameters
        ----------
        self

        Returns
        -------
        None

        """
        self.__is_hotspot = True if self.__gene in IMPACTV6_HOTSPOTS and self.__is_coding else False

    def __set_is_kinase(self):
        r"""
        Sets __is_kinase flag based on if gene is in IMAPCTV6_KINASE_TARGETS in constants file

        Parameters
        ----------
        self

        Returns
        -------
        None

        """
        self.__is_kinase = True if self.__gene in IMAPCTV6_KINASE_TARGETS and self.__is_coding else False

    def __set_tumor_suppressor(self):
        r"""
        Sets __is_tumour_suppressor flag based on if gene is in IMPACTV6_TUMOR_SUPPRESSOR in constants file
        and gene is in coding

        Parameters
        ----------
        self

        Returns
        -------
        None

        """
        self.__is_tumour_suppressor = True if self.__gene in IMPACTV6_TUMOR_SUPPRESSOR and self.__is_coding else False

    def __replace_ensembl_tx(self):
        r"""
        Genes RECQL4, PRIM2, IKZF1 have RefSec transcript mapping issue hence we use Ensembl transcript
        This method replaces Ensembl with RefSec transcript for these genes

        Parameters
        ----------
        self

        Returns
        -------
        tx : string
            RefSec transcript

        """
        __tx = self.__transcript
        if not __tx.startswith("ENST"):
            return __tx
        __lk_up = {"ENST00000331340.8":"NM_006060.6", "ENST00000428558.2":"NM_004260.3", "ENST00000607273.1":"NM_000947.5"} 
        return __lk_up.get(__tx,"")

    def __set_transcript(self):
        r"""
        Read transctipt for the gene from ref_flat_canonical in constants file

        Parameters
        ----------
        self

        Returns
        -------
        None

        """
        try:
            self.__transcript = ref_flat_canonical[self.__gene]
        except KeyError as e:
            self.__logger.error("Cannot find canonical transcript for {}".format(self.__gene))
            self.__logger.error(e)

        try:
            self.__transcript, self.__cdna = get_cdna_pos(self,self.__logger)
            self.__transcript = self.__replace_ensembl_tx()
            self.__transcript = self.__transcript.split(".")[0]
        except Exception as e:
            self.__logger.warning(e)

    def __expand(self):
        r"""
        This method is used to form different variables of breakpoint class

        Parameters
        ----------
        self

        Returns
        -------
        None

        """
        try:
            self.__set_transcript()
            self.__set_strand()
            self.__set_is_coding()
            self.__set_is_panel()
            self.__set_is_hotspot()
            self.__set_is_kinase()
            self.__set_tumor_suppressor()
        except Exception as e:
            self.__logger.error("Exception in expand method")
            self.__logger.error(e)
        
    def set_is_entire_kinase(self, val): 
        r"""
        Public method to __is_tumour_suppressor

        Parameters
        ----------
        self
        val : int
            0 : entire kinase not invloved
            1 : entire kinase invloved
            -1 : entire kinase partial

        Returns
        -------
        None

        """
        self.__is_entire_kinase = val

    def set_first_exon(self,val): 
        r"""
        Public method to set __first_exon

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__first_exon = int(val)

    def set_last_exon(self,val):  
        r"""
        Public method to set __last_exon

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__last_exon = int(val)

    def set_start_pos(self,val):  
        r"""
        Public method to set __start_pos

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__start_pos = int(val)

    def set_stop_pos(self,val):  
        r"""
        Public method to set __stop_pos

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__stop_pos = int(val)

    def set_exon(self,val):  
        r"""
        Public method to set __exon

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__exon = int(val)

    def set_intron(self,val):  
        r"""
        Public method to set __intron

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__intron = val

    def set_site(self,val):  
        r"""
        Public method to set __site

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__site = val

    def set_variant_site_1(self,val):  
        r"""
        Public method to set __variant_site_1

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__variant_site_1 = val

    def set_variant_site_2(self,val):  
        r"""
        Public method to set __variant_site_2

        Parameters
        ----------
        self
        val : string

        Returns
        -------
        None

        """
        self.__variant_site_2 = val
        
    def get_chrom(self): 
        r"""
        Public method to get __chrom

        Parameters
        ----------
        self

        Returns
        -------
        __chrom : string

        """
        return self.__chrom

    def get_pos(self): 
        r"""
        Public method to get __pos

        Parameters
        ----------
        self

        Returns
        -------
        __pos : int

        """ 
        return self.__pos

    def get_gene(self): 
        r"""
        Public method to get __gene

        Parameters
        ----------
        self

        Returns
        -------
        __gene : string

        """
        return self.__gene

    def get_desc(self):  
        r"""
        Public method to get __desc

        Parameters
        ----------
        self

        Returns
        -------
        __desc : string

        """
        return self.__desc

    def get_cdna(self):  
        r"""
        Public method to get __chrom

        Parameters
        ----------
        self

        Returns
        -------
        __chrom : string

        """
        return self.__cdna

    def get_transcript(self):  
        r"""
        Public method to get __transcript

        Parameters
        ----------
        self

        Returns
        -------
        __transcript : string

        """
        return self.__transcript

    def get_is_coding(self):  
        r"""
        Public method to get __chrom

        Parameters
        ----------
        self

        Returns
        -------
        __chrom : string

        """
        return self.__is_coding

    def get_is_hotspot(self):  
        r"""
        Public method to get __is_hotspot

        Parameters
        ----------
        self

        Returns
        -------
        __is_hotspot : bool

        """
        return self.__is_hotspot

    def get_is_kinase(self):  
        r"""
        Public method to get __is_kinase

        Parameters
        ----------
        self

        Returns
        -------
        __is_kinase : bool

        """
        return self.__is_kinase

    def get_is_kinase_hotspot(self):  
        r"""
        Public method to get __is_hotspot

        Parameters
        ----------
        self

        Returns
        -------
        __is_hotspot : bool

        """
        return self.__is_kinase and self.__is_hotspot

    def get_is_panel(self):  
        r"""
        Public method to get __is_panel

        Parameters
        ----------
        self

        Returns
        -------
        __is_panel : bool

        """
        return self.__is_panel

    def get_is_panel_and_coding(self):  
        r"""
        Public method to get __is_coding and __is_panel

        Parameters
        ----------
        self

        Returns
        -------
        __is_panel_coding : bool

        """
        return self.__is_panel and self.__is_coding

    def get_tumor_suppressor(self):  
        r"""
        Public method to get __is_tumour_suppressor

        Parameters
        ----------
        self

        Returns
        -------
        __is_tumour_suppressor : bool

        """
        return self.__is_tumour_suppressor

    def get_is_entire_kinase(self):  
        r"""
        Public method to get __is_entire_kinase

        Parameters
        ----------
        self

        Returns
        -------
        __is_entire_kinase : int
            0 : entire kinase not invloved
            1 : entire kinase invloved
            -1 : entire kinase partial

        """
        return self.__is_entire_kinase

    def get_first_exon(self):  
        r"""
        Public method to get __first_exon

        Parameters
        ----------
        self

        Returns
        -------
        __first_exon : int

        """
        return self.__first_exon

    def get_last_exon(self):   
        r"""
        Public method to get __last_exon

        Parameters
        ----------
        self

        Returns
        -------
        __last_exon : int

        """
        return self.__last_exon
        
    def get_start_pos(self):   
        r"""
        Public method to get __start_pos

        Parameters
        ----------
        self

        Returns
        -------
        __start_pos : int

        """
        return self.__start_pos

    def get_stop_pos(self):   
        r"""
        Public method to get __stop_pos

        Parameters
        ----------
        self

        Returns
        -------
        __stop_pos : int

        """
        return self.__stop_pos

    def get_exon(self):   
        r"""
        Public method to get __exon

        Parameters
        ----------
        self

        Returns
        -------
        __exon : int

        """
        return self.__exon

    def get_intron(self):   
        r"""
        Public method to get __intron

        Parameters
        ----------
        self

        Returns
        -------
        __intron : string

        """
        return self.__intron

    def get_site(self):   
        r"""
        Public method to get __site

        Parameters
        ----------
        self

        Returns
        -------
        __site : string

        """
        return self.__site

    def get_strand(self):   
        r"""
        Public method to get __strand

        Parameters
        ----------
        self

        Returns
        -------
        __strand : string

        """
        return self.__strand

    def get_variant_site_1(self):   
        r"""
        Public method to get __variant_site_1

        Parameters
        ----------
        self

        Returns
        -------
        __variant_site_1 : string

        """
        return self.__variant_site_1

    def get_variant_site_2(self):   
        r"""
        Public method to get __variant_site_2

        Parameters
        ----------
        self

        Returns
        -------
        __variant_site_2 : string

        """
        return self.__variant_site_2

    
    def check_is_entire_kinase_involved(self):    
        r"""
        Public method to check if entire kinase is involved

        Parameters
        ----------
        self

        Returns
        -------
        _is_entire_kinase_involved : bool

        Note:
        ---------
        Returns true if __is_entire_kinase = 1

        """
        return int(self.__is_entire_kinase)  == 1

    def check_is_kinase_partial(self):    
        r"""
        Public method to check if entire kinase is partial

        Parameters
        ----------
        self

        Returns
        -------
        _is_entire_kinase_partial : bool

        Note:
        ---------
        Returns true if __is_entire_kinase = -1

        """
        return int(self.__is_entire_kinase)  == -1

    def check_is_kinase_not_involved(self):    
        r"""
        Public method to check if entire kinase is not involved

        Parameters
        ----------
        self

        Returns
        -------
        _is_entire_kinase_not_involved : bool

        Note:
        ---------
        Returns true if __is_entire_kinase = 0

        """
        return int(self.__is_entire_kinase)  == 0
