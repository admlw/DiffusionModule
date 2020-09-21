#INPUT_FILE="/pnfs/uboone/scratch/users/amogan/v08_00_00_25/run3_g1_crt_throughgoing_SCE/diffusion_filter/33654394_39/PhysicsRun-2017_12_11_14_39_29-0014224-00056_20180811T061028_ext_bnb_2_20181213T074227_optfilter_20181229T033215_reco1_postwcct_postdl_20181229T034029_d45fab9a-3b65-4b8d-b8ea-7f766d521fb9.root"
#OUTPUT_FILE=mp.root

#WIRE_PRODUCER="butcher"
#RUN=14146
#MIN_EVENT=1477
#MAX_EVENT=1479
#MIN_CHANNEL=6900
#MAX_CHANNEL=7000

INPUT_FILE="/uboone/app/users/amogan/diffusion_mcc9/workdir/prod_muminus_fwdgoing_custom_20200831T160954_gen_20200831T161027_g4_20200831T161921_detsim.root"
OUTPUT_FILE="out.root"

#RAWD_PRODUCER="wcNoiseFilter"
RAWD_PRODUCER="driftWC:orig:Detsim"
RUN=1
MIN_EVENT=1
MAX_EVENT=1
NTICKS=6400
MIN_CHANNEL=5500
MAX_CHANNEL=5500

#if [ -z "$4" ]
#then
#    MAX_CHANNEL=$3
#fi

make all
./GetWire "$INPUT_FILE" "$OUTPUT_FILE" "$WIRE_PRODUCER" "$RUN" "$MIN_EVENT" "$MAX_EVENT" "$MIN_CHANNEL" "$MAX_CHANNEL"

rm GetWire
