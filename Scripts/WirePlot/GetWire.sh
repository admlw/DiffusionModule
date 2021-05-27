#INPUT_FILE="/pnfs/uboone/scratch/users/amogan/v08_00_00_25/run3_g1_crt_throughgoing_SCE/diffusion_filter/33654394_39/PhysicsRun-2017_12_11_14_39_29-0014224-00056_20180811T061028_ext_bnb_2_20181213T074227_optfilter_20181229T033215_reco1_postwcct_postdl_20181229T034029_d45fab9a-3b65-4b8d-b8ea-7f766d521fb9.root"
#OUTPUT_FILE=mp.root

#WIRE_PRODUCER="butcher"
#RUN=14146
#MIN_EVENT=1477
#MAX_EVENT=1479
#MIN_CHANNEL=6900
#MAX_CHANNEL=7000

#INPUT_FILE="/uboone/app/users/amogan/diffusion_mcc9/workdir/prod_muminus_fwdgoing_custom_20200831T160954_gen_20200831T161027_g4_20200831T161921_detsim.root"
#INPUT_FILE="/pnfs/uboone/persistent/users/alister1/diffusion_crt_data_forPaper/filtered/PhysicsRun-2018_3_10_5_1_28-0015425-00174_20180817T100002_ext_bnb_4_20181216T215102_optfilter_20181230T152648_reco1_postwcct_postdl_20181230T155854_20_d1639faf-e9c2-4f9c-a371-1d3a7a480bb3.root"
INPUT_FILE="/pnfs/uboone/persistent/diffusion_analysis_final_files/filtered_data_files/PhysicsRun-2018_3_10_5_1_28-0015425-00174_20180817T100002_ext_bnb_4_20181216T215102_optfilter_20181230T152648_reco1_postwcct_postdl_20181230T155854_20_d1639faf-e9c2-4f9c-a371-1d3a7a480bb3.root"
#INPUT_FILE="/pnfs/uboone/persistent/diffusion_analysis_final_files/filtered_data_files/PhysicsRun-2018_6_6_12_30_15-0017013-00046_20181001T045313_ext_bnb_1_20181224T172511_optfilter_20200519T145645_reco1_postwcct_postdl_20200519T151009_2_eafdaf6b-29cc-42ef-9dd0-12d00f8841b1.root"
OUTPUT_FILE="out.root"

#RAWD_PRODUCER="wcNoiseFilter"
WIRE_PRODUCER="butcher"
RUN=14379
MIN_EVENT=66
MAX_EVENT=66
NTICKS=6400
MIN_CHANNEL=7000
MAX_CHANNEL=7001

#if [ -z "$4" ]
#then
#    MAX_CHANNEL=$3
#fi

make all
./GetWire "$INPUT_FILE" "$OUTPUT_FILE" "$WIRE_PRODUCER" "$RUN" "$MIN_EVENT" "$MAX_EVENT" "$MIN_CHANNEL" "$MAX_CHANNEL"

rm GetWire
