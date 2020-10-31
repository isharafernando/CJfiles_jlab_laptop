#/bin/bash

 ###### This is important for mac terminal: ###############
 [ -f ~/.profile ] && source ~/.profile || echo " "

 ###### Recompile Cooker: ##################################
 cd build
 sudo make install
 wait
 cd ..

echo "here's that gem control..."
# midas to root conversion

#cooker recipes/midasconverter.xml /home/ishara/MUSE/midfiles/run0$1.mid /home/ishara/MUSE/mid2rootfiles/mid2root0$1.root

# GEM cluster finding:
#cooker recipes/gemControl.xml  $ROOTFILES/run$1.root  $ROOTFILES/cookedGEMrun$1.root
#cooker recipes/midasconverter.xml /home/ishara/MUSE/GEM-Data/run0$1.mid /home/ishara/MUSE/cookedfiles/cookedGEMrun0$1.root

#cooker recipes/gemControl.xml /home/ishara/MUSE/mid2rootfiles/mid2root0$1.root /home/ishara/MUSE/cookedfiles/cookedGEMrun0$1.root

# GEM tracker:
#cooker recipes/MUSEteleTracker_gemeff.xml $ROOTFILES/cookedGEMrun$1.root $ROOTFILES/GEM_tracked$1.root

#cooker recipes/MUSEteleTracker_gemeff.xml /home/ishara/MUSE/cookedfiles/cookedGEMrun0$1.root /home/ishara/MUSE/trackedfiles/GEM_tracked0$1.root

#1 midas to root conversion
#cooker recipes/midasconverter.xml /home/ishara/MUSE/midfiles/run0$1.mid /home/ishara/MUSE/mid2rootfiles/mid2root0$1.root

#2 GEM cluster finding:
#cooker recipes/gemControl.xml /home/ishara/MUSE/mid2rootfiles/mid2root0$1.root /home/ishara/MUSE/cookedfiles/cookedGEMrun0$1.root

#3 BH rootfile preparation
#cooker recipes/BH.xml /home/ishara/MUSE/mid2rootfiles/mid2root0$1.root /home/ishara/MUSE/BHcookedfiles/cooked_BH0$1.root

#4 GEM tracker:
cooker recipes/MUSEteleTracker_gemeff.xml /home/ishara/MUSE/BHcookedfiles/cooked_BH0$1.root:/home/ishara/MUSE/cookedfiles/cookedGEMrun0$1.root /home/ishara/MUSE/trackedfiles/GEM_BH_tracked0$1.root
