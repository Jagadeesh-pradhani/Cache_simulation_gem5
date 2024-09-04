--------------------------------------------------------------------------------
#Custom protocol Custom_MESI_Two_Level

#Create a new folder in gem5/src 

mkdir ~/gem5/src/REU

#Create new sm and slicc files in that folder.
Custom_MESI_Two_Level.slicc
Custom_MESI_Two_Level-msg.sm
Custom_MESI_Two_Level-L1cache.sm
Custom_MESI_Two_Level-L2cache.sm
Custom_MESI_Two_Level-dir.sm
Custom_MESI_Two_Level-dma.sm

#Building
scons build/X86_Custom_MESI_Two_Level/gem5.opt --default=X86 PROTOCOL=Custom_MESI_Two_Level