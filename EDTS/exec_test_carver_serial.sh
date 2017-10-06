#!/bin/bash

#set -x
#===============================================================================
#                    User control variables
#===============================================================================

# This is the unique identifier tag for this test environment
# It should indicate the revision number branched from, the users initials
# and a test number.  For instance if you branched from r81, your initials
# are rk and you tried 4 iterations of this test before your commits
# were verified, you would use:

VERSION="r85v2ghub"

# FILE PATHS TO YOUR THREE EXECUTABLES

#MAIN_EXE_PATH=$VSC_SCRATCH_VO/manfredo/EDTSmine/binaries/ed_2.1-opt_raichu_MAIN
#TEST_EXE_PATH=$VSC_SCRATCH_VO/manfredo/EDTSmine/binaries/ed_2.1-opt_raichu_TEST
#DBUG_EXE_PATH=$VSC_SCRATCH_VO/manfredo/EDTSmine/binaries/ed_2.1-dbg_raichu_DBUG
MAIN_EXE_PATH=/user/scratch/gent/vsc416/vsc41690/ED2_master/build/ed_2.1-opt_golett
TEST_EXE_PATH=/user/scratch/gent/vsc416/vsc41690/ED_redux/build/ed_2.1-opt_golett
DBUG_EXE_PATH=/user/scratch/gent/vsc416/vsc41690/ED_redux/build/ed_2.1-dbg_golett

# Provide the path where the test-suit driver archive is stored

DATAPATH="$VSC_DATA_VO/manfredo/edts_datasets"

# Decide on a "rapid" (2 year) or "long" (300 year)
# THE long tests will take a really really long time, beware

TESTTYPE="rapid"

# Decide on which sites to use.  You should ultimately
# have results for all of them, but this is helpfull if some
# sites give you problems after the first try
# Y for Yes and N for No

USE_M34=Y     #Manaus km 34 SOI
USE_S67=Y     #Santarem km 67 SOI
USE_HAR=Y     #Harvard forest SOI
USE_PDG=Y     #Pe de Gigante SOI 
USE_TON=N     #Tonzi SOI
USE_CAX=Y     #Caxuana SOI
USE_TNF=Y     #Tapajos National Forest SOI
USE_ATA=Y     #Atacama Desert SOI
USE_GYF=Y     #Paracou
USE_PET=Y     #Petrolina SOI
USE_HIP=Y     #Petrolina High Frequency Detailed Short SOI
USE_HIM=Y     #Manaus High Frequency Detailed Short SOI
USE_RJG=N     #Gridded 12x12 simulation centerd on Reserva Jaru

# How many cores do you want to use for the gridded simulations
# Currently there are 3 (RJG-MAIN RJG-TEST and RJG-DBUG)

NNODES=3
PPN=8

# Give an explanation of the tests.  Explain what the commits had involved.

TEST_DESCRIPTION="This is a test of the LIANAS + CLEANING pull request #200. Conducted by Manfredo"

# The identifier may had indicated which version you branched from, but indicate it here
# also

VERSION_BRANCHED_FROM='r84'

# Who is running this test?
TESTER_NAME='Manfredo'

# Who was the developer(s) that actually made the changes to the code that is being tested?
COMMITTER_NAME='Manfredo'



#===============================================================================


echo ""
echo ""
echo "========================================================================="
echo "                 Starting the EDM Dev Test Suit  (EDTS)                  "
echo "========================================================================="
echo ""

# DEFINE SOME RUNTIME VARIABLES FOR SOI's

declare -a USE_SITE=( $USE_M34 $USE_S67 $USE_HAR $USE_PDG $USE_TON $USE_CAX \
    $USE_TNF $USE_ATA $USE_PET $USE_GYF $USE_HIP $USE_HIM $USE_RJG)
declare -a SITEID=(m34 s67 har pdg ton cax tnf ata pet gyf hip him rjg)
declare -a SITEPFX=(M34 S67 HAR PDG TON CAX TNF ATA PET GYF HIP HIM RJG)

# SOI DEBUG 
declare -a IYEARAS=(1500 1500 2007 1500 2000 2000 2002 1500 2005 2005 1500 1500 2008)
declare -a IYEARZS=(1502 1502 2009 1502 2002 2002 2004 1502 2007 2007 1500 1500 2008)
declare -a INITMDS=(5    0    6    0    6    0    5    0    0    0    6    5    5)

#TYPE OF RUN
declare -a RUNTYPE=(reg reg reg reg reg reg reg reg reg reg hif hif grd)

# HI FREQUENCY DETAILED RUNS (HI-PET and HI-M34)
declare -a IDATEAH=(21 01)
declare -a IDATEZH=(28 08)

UNITSTATE=2

let NPROC=NNODES*PPN


for i in ${!SITEID[@]}
do
    if [ ${USE_SITE[i]} == "Y" ]; then
        echo "PROCESSING ${RUNTYPE[i]}: "${SITEID[i]}
    elif [ ${USE_SITE[i]} == "N" ]; then
        echo "SKIPPING ${RUNTYPE[i]}: "${SITEID[i]}
    else
        echo "IMPROPER USE SPECIFIER: "${SITEID[i]}
        echo "STOPPING"
        #exit
    fi
done


#===============================================================================
#===============================================================================


# Update version string to contain the test type and write some useful arrays
VERSION=${VERSION}${TESTTYPE}
EXE=( "./ed_2.1-main" "./ed_2.1-test" "./ed_2.1-dbug" )
RUN_KIND=( MAIN TEST DBUG )
EXE_PATH=( $MAIN_EXE_PATH $TEST_EXE_PATH $DBUG_EXE_PATH )

# Create the working folder

echo ""
echo ""
echo "============================================="
echo "  GENERATE PRISTINE DIRECTORY: ${VERSION} ???"
echo ""
echo " ALL DATA WILL BE LOST!"
echo " ANSWER [Y/N]"
echo "============================================="


read flusher

if [ ${flusher} == "Y" ]; then
    echo "FLUSHING"
    sleep 2
    rm -rf ${VERSION}
    mkdir -p ${VERSION}

    # Generate the symbolic link to edts_datasets
    for  kind in ${!RUN_KIND[@]}; do
        ln -s ${EXE_PATH[${kind}]} ${VERSION}/${EXE[${kind}]}
    done

    ln -s ${DATAPATH} ${VERSION}/edts_datasets
    cp Templates/*xml ${VERSION}/

else
    if [ ${flusher} == "N" ]; then
        echo "WILL NOT FLUSH"
        sleep 2
    else
        echo "NEITHER Y OR N, EXITING"
        exit
    fi
fi


# Write out the text info to the report directory

echo '<?xml version="1.0" ?>'                                       > ${VERSION}/test_text.xml
echo '<description> '                                              >> ${VERSION}/test_text.xml
echo '<branch_version> '$VERSION_BRANCHED_FROM' </branch_version>' >> ${VERSION}/test_text.xml
echo '<tester_name> '$TESTER_NAME' </tester_name>'                 >> ${VERSION}/test_text.xml
echo '<committer_name> '$COMMITTER_NAME' </committer_name>'        >> ${VERSION}/test_text.xml
echo '<test_description> '$TEST_DESCRIPTION' </test_description>'  >> ${VERSION}/test_text.xml
echo '</description>'                                              >> ${VERSION}/test_text.xml


# Loop over the different SOI CASES
# ==============================================================================

ntasks=0
nnodes=1

for i in ${!SITEID[@]}
do
    if [ ${USE_SITE[i]} == "Y" ]; then
        let ntasks=ntasks+3
    fi
done

if [ ${ntasks} > 0 ]; then
    echo "TOTAL NUMBER OF SOI TASKS = "${ntasks}
    let nnodes=nnodes+`awk  'BEGIN { rounded = sprintf("%.0f", '${ntasks}/24' ); print rounded }'`
    echo "TOTAL NUMBER OF NODES = "${nnodes}
    let mppwidth=`awk  'BEGIN { rounded = sprintf("%.0f", '${nnodes}*24' ); print rounded }'`
    echo "MPPWIDTH = "${mppwidth}

else
    echo "NO SOI TASKS"
fi




echo '#!/bin/sh' > ${VERSION}/submit_batch.sh
chmod +x ${VERSION}/submit_batch.sh

hif_count=0
for i in ${!SITEID[@]}; do

    if [ ${USE_SITE[i]} == "Y" ]; then

        TEMPLATE=Templates/ED2IN-${SITEPFX[i]}-MAIN
        echo ""
        echo "Processing "${SITEID[i]}
        echo ""

        for j in ${!RUN_KIND[@]}; do

            kind=${RUN_KIND[j]}
            skind=$(echo "$kind" | awk '{print tolower($0)}')
            FILE=${VERSION}/ED2IN-${SITEPFX[i]}-$kind

            echo "Control File: "$FILE

            cp $TEMPLATE $FILE

            # Modify the ED2IN RUNTYPE
            # ========================
            sed -i '/NL%RUNTYPE/c\   NL%RUNTYPE = '\'INITIAL\' $FILE

            # Modify the ED2IN IED_INIT_MODE
            # ==============================
            sed -i '/NL%IED_INIT_MODE/c\   NL%IED_INIT_MODE = '${INITMDS[i]} $FILE

            # Modify the ED2IN start years
            # ============================
            sed -i '/NL%IYEARA/c\   NL%IYEARA   = '${IYEARAS[i]} $FILE

            # Modify the ED2IN end years
            # ==========================
            sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARZS[i]} $FILE

            # Modify the ED2IN files to point to the desired output directories   
            # =================================================================   
            sed -i '/NL%FFILOUT/c\   NL%FFILOUT = '\''F_'$skind'_'${SITEID[i]}'/'$skind'_'${SITEID[i]}''\' $FILE
            sed -i '/NL%SFILOUT/c\   NL%SFILOUT = '\''S_'$skind'_'${SITEID[i]}'/'$skind'_'${SITEID[i]}''\' $FILE

            # Modify the state file output frequency
            # ======================================
            sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = '$UNITSTATE $FILE


            if [ ${RUNTYPE[i]} != "reg" ]; then

                # Modify the ED2IN end years
                # ==========================
                sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '${IYEARAS[i]} $FILE

                if [ ${RUNTYPE[i]} == "hif" ]; then

                    # Modify the ED2IN start dates
                    # ============================
                    sed -i '/NL%IDATEA/c\   NL%IDATEA   = '${IDATEAH[$hif_count]} $FILE

                    # Modify the ED2IN end date
                    # ==========================
                    sed -i '/NL%IDATEZ/c\   NL%IDATEZ   = '${IDATEZH[$hif_count]} $FILE

                    # Modify the state file output frequency
                    # ======================================
                    sed -i '/NL%UNITSTATE/c\   NL%UNITSTATE  = 3' $FILE

                    if [ $kind == ${RUN_KIND[-1]} ]; then
                        ((hif_count++))
                    fi 

                fi

            fi

            if [ $kind == "DBUG" ] && [ ${RUNTYPE[i]} == "reg" ]; then

                # Modify the ED2IN end years
                # ==========================
                sed -i '/NL%IYEARZ/c\   NL%IYEARZ   = '$((IYEARAS[$i]+1)) $FILE

            fi

            # Manfredo: I removed the IOPTINT option, here we remove it from the ED2IN  
            if [ $kind != "MAIN" ]; then

                sed -i '/NL%IOPTINPT/d' $FILE
                sed -i '/NL%IBIGLEAF/d' $FILE

            fi

            # Reset and flush the output folders
            # ====================================
            mkdir -p ${VERSION}/F_${skind}_${SITEID[i]}
            mkdir -p ${VERSION}/S_${skind}_${SITEID[i]}

            # CREATE THE PBS FILE
            # ===================
            PBSFILE=${VERSION}/batch_${skind}_${SITEID[i]}
            echo "#PBS -l walltime=18:00:00"                          > $PBSFILE
            if [ ${RUNTYPE[i]} == "grd" ]; then
                echo -e "#PBS -l nodes=${NNODES}:ppn=${PPN}"           >> $PBSFILE
            fi
            echo "#PBS -l pvmem=3GB"                                 >> $PBSFILE
            echo '#PBS -N '${kind}'_'${SITEID[i]}                    >> $PBSFILE
            echo '#PBS -e '${skind}'_'${SITEID[i]}'.$PBS_JOBID.err'  >> $PBSFILE
            echo '#PBS -o '${skind}'_'${SITEID[i]}'.$PBS_JOBID.out'  >> $PBSFILE
            echo "#PBS -V"                                           >> $PBSFILE
            echo ""                                                  >> $PBSFILE
            echo "module load HDF5; ulimit -s unlimited"             >> $PBSFILE
            echo ""                                                  >> $PBSFILE
            echo 'cd $PBS_O_WORKDIR'                                 >> $PBSFILE
            echo ""                                                  >> $PBSFILE
            if [ ${RUNTYPE[i]} == "grd" ]; then
                echo -e "mpirun -n ${NPROC} \c"                        >> $PBSFILE
            fi
            echo "${EXE[$j]} -f ED2IN-${SITEPFX[i]}-${kind}"         >> $PBSFILE
            echo ""                                >> ${VERSION}/submit_batch.sh
            echo "qsub batch_${skind}_${SITEID[i]}">> ${VERSION}/submit_batch.sh

        done


        echo "IO/Err Files: "${VERSION}/${SITEID[i]}"*"
        echo ""

    elif [ ${USE_SITE[i]}=="N" ]; then

        if [ ${RUNTYPE[i]} == "hif" ]; then
            ((hif_count++))
        fi
        echo "Skipping"

    else

        echo "YOUR USE_SITE VECTOR HAS AN ERROR"
        exit

    fi

done
