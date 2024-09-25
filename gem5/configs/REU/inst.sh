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



build/X86/gem5.opt --debug-flags=ProtocolTrace --debug-file=output1.out configs/deprecated/example/se.py --cmd="/home/jagadeesh/gem5/SPEC.Small/mcf/mcf;/home/jagadeesh/gem5/SPEC.Small/dct/dct"  --options="" -n 2 --num-dirs=2 --ruby --network=garnet --topology=Mesh_XY --mesh-rows=1 --link-latency=1 --router-latency=1 --cpu-type=O3CPU --l1d_size=32kB --l1i_size=32kB --l2_size=4MB --mem-size=512MB --l1d_assoc=2 --l1i_assoc=2 --l2_assoc=8