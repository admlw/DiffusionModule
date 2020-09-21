#INPUT_FILE="/uboone/app/users/amogan/diffusion_mcc9/workdir/prod_muminus_fwdgoing_custom_20200831T160954_gen_20200831T161027_g4_20200902T170424_detsim_DTwaywayup.root"
#INPUT_FILE="/uboone/app/users/amogan/diffusion_mcc9/workdir/prod_muminus_fwdgoing_custom_20200831T160954_gen_20200831T161027_g4_20200901T204839_detsim.root"
#INPUT_FILE="/uboone/app/users/amogan/diffusion_mcc9/workdir/prod_muminus_fwdgoing_custom_20200903T202407_gen_20200903T202449_g4_20200903T204103_detsim_DToff.root"
INPUT_FILE="/uboone/app/users/amogan/diffusion_mcc9/workdir/prod_muminus_fwdgoing_custom_20200903T202407_gen_20200903T202449_g4_20200903T212517_detsim_DTwaywayup.root"
OUTPUT_FILE="/uboone/data/users/amogan/wireplot_outfiles_diffusion/out_DT_way_way_up_thetaXZ_20.root"

#RAWD_PRODUCER="wcNoiseFilter"
RAWD_PRODUCER="driftWC:orig:DetsimDTwaywayup"
#RAWD_PRODUCER="driftWC:orig:DetsimDToff"
RUN=1
MIN_EVENT=0
MAX_EVENT=10
NTICKS=6400
#MIN_CHANNEL=$3
#MAX_CHANNEL=$4
MIN_CHANNEL=0
MAX_CHANNEL=8256

#if [ -z "$4" ]
#then
#    MAX_CHANNEL=$3
#fi

make all
./GetWaveform "$INPUT_FILE" "$OUTPUT_FILE" "$RAWD_PRODUCER" "$RUN" "$MIN_EVENT" "$MAX_EVENT" "$MIN_CHANNEL" "$MAX_CHANNEL" "$NTICKS"

rm GetWaveform
