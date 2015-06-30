#!/usr/bin/env python
###############################################################################
#                                                                             #
#    viewInterface.py                                                         #
#                                                                             #
#    Interface for accessing hit data                                         #
#                                                                             #
#    Copyright (C) Joshua Daly                                                #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Joshua Daly"
__copyright__ = "Copyright 2015"
__credits__ = ["Joshua Daly"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Joshua Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# system includes
import sys

# local includes
import fileparser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ViewInterface(Interface):
    """Use this interface for querying the TrackM DB"""

    def __init__(self,
                 hit_data, 
                 transfer_group_file,
                 contam_pidsqids_file
                 ):
        self.TG                             = TFP.GroupData(transfer_group_file)
        self.HD                             = TFP.HitData(hit_data)
        self.clean_analysis                 = contam_pidsqids_file
        if self.clean_analysis:
            self.CP                         = TFP.ContaminatedPidsqids(contam_pidsqids_file)
        


###############################################################################
###############################################################################
###############################################################################
###############################################################################

