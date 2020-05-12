# Takes artroot file as input, along with run, event, and channel numbers
# Currently using /pnfs/uboone/scratch/users/amogan/v08_00_00_25/run3_g1_crt_filterVol/diffusion_filter/7517332_52/PhysicsRun-2018_3_21_12_6_36-0015670-00002_20190124T070718_ext_bnb_2_20190130T230330_optfilter_20190217T061552_reco1_postwcct_20190217T061734_20191028_fb52f80f-ff45-419f-992a-197568a886e6.root
# Run: 14436 subrun: 99 event: 4973 track_start_z: 483.89285 

INPUT_FILE="/pnfs/uboone/scratch/users/amogan/v08_00_00_25/run3_g1_crt_filterVol/diffusion_filter/7517332_52/PhysicsRun-2018_3_21_12_6_36-0015670-00002_20190124T070718_ext_bnb_2_20190130T230330_optfilter_20190217T061552_reco1_postwcct_20190217T061734_20191028_fb52f80f-ff45-419f-992a-197568a886e6.root"
OUTPUT_FILE=out.root

#WIRE_PRODUCER="caldata"
WIRE_PRODUCER="butcher"
RUN=14436
MIN_EVENT=4972
MAX_EVENT=4974
MIN_CHANNEL=6150
MAX_CHANNEL=8000

if [ -z "$4" ]
then
    MAX_CHANNEL=$3
fi

make all
./GetWire "$INPUT_FILE" "$OUTPUT_FILE" "$WIRE_PRODUCER" "$RUN" "$MIN_EVENT" "$MAX_EVENT" "$MIN_CHANNEL" "$MAX_CHANNEL"

rm GetWire
