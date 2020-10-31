/*
 * Definitions of the gemControl class and methods
 * Noah Wuerfel nwuerfel@umich.edu 10-12-18 goblu
 * ~AP AP AP AP~
 */

#include <gemControl.h>
#include <iostream>
#include <fstream>
#include <cmath>

gemControl::gemControl(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
    // gem's labeled by ID from upstream to downstream
    // cleaned up in destructor
    for(std::vector<gem*>::size_type i=0; i < NUM_GEMS*GEM_BANK_NUM; i++){
        gem* gem_to_add = new gem(i);
        if(!gem_to_add){
            std::cout << "failed to create gem with id: " << i
                << std::endl;
            exit(-1);
        }
        gem_list.push_back(gem_to_add);
    }
}

// clean up
gemControl::~gemControl(){
    for (int i = 0; i<GEM_BANK_NUM; i++){
      std::cout << "cleaning up..."<<GEM_bank[i]<<"\n";
      std::cout << "cleaning up blob reader...\n";
      delete binBlobReader[i];
    }
    std::cout << "cleaning up gems...\n";
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        assert(gem_list[i]!=NULL);
        delete gem_list[i];
    }
    delete gemo;

}

Long_t gemControl::defineHistograms(){
    std::cout << "here it is, the big gemControl...\n";
    std::vector<TH2D*> temp_apv_hist_vec;
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        const char* name = gem_list[i]->name;

        TH2D* gem_h2_x = dH2(Form("%s/adcSpec_gem_%i_x",name,(int)i),
            Form("%s GEM X spectrum ;APV Strip Number; ADC Value",name), NUM_CHAN_PER_AXIS,
            SPEC_RANGE_MIN, NUM_CHAN_PER_AXIS,NUM_CHAN_PER_AXIS,
            RANGE_MIN, RANGE_MAX);
        TH2D* gem_h2_y = dH2(Form("%s/adcSpec_gem_%i_y",name,(int)i),
            Form("%s GEM Y spectrum; APV Strip Number; ADC Value",name), NUM_CHAN_PER_AXIS,
            SPEC_RANGE_MIN, NUM_CHAN_PER_AXIS,NUM_CHAN_PER_AXIS,
            RANGE_MIN, RANGE_MAX);

        TH2D* gem_h2_x_cmode = dH2(Form("%s/adcSpec_gem_%i_x_cmode",name,(int)i),
            Form("%s GEM X spectrum cmode; APV Strip Number; ADC Value",name), NUM_CHAN_PER_AXIS,
            SPEC_RANGE_MIN, NUM_CHAN_PER_AXIS,NUM_CHAN_PER_AXIS,
            RANGE_MIN, RANGE_MAX);
        TH2D* gem_h2_y_cmode = dH2(Form("%s/adcSpec_gem_%i_y_cmode",name,(int)i),
            Form("%s GEM Y spectrum cmode;APV Strip Number; ADC Value",name), NUM_CHAN_PER_AXIS,
            SPEC_RANGE_MIN, NUM_CHAN_PER_AXIS,NUM_CHAN_PER_AXIS,
            RANGE_MIN, RANGE_MAX);

        TH2D* gem_hitmap = dH2(Form("%s/hitmap_gem_%i",name,(int)i),
            Form("%s GEM Hitmap; X Position(APV Strip Number); Y Position(APV Strip Number)",name), NUM_CHAN_PER_AXIS,
            SPEC_RANGE_MIN,NUM_CHAN_PER_AXIS, NUM_CHAN_PER_AXIS,
            SPEC_RANGE_MIN, NUM_CHAN_PER_AXIS);
        TH2D* gem_clustermap = dH2(Form("%s/clustermap_gem_%i",name,(int)i),
            Form("%s GEM Clustermap; X Position(mm); Y Position(mm)",name), 100, -50, 50, 100, -50, 50);

       TH1D* gem_cluster_mult = dH1(Form("%s/cluster_multiplicity_gem_%i",name,(int)i),
            Form("%s Cluster Multiplicity",name), 100, -0.5, 99.5);

        for(std::vector<apv*>::size_type j=0;
            j<gem_list[i]->apv_list.size(); j++){
            uint32_t apv_id = gem_list[i]->apv_list[j]->id;
            uint32_t num_skip = gem_list[i]->apv_list[j]->\
                num_chan_to_skip;
            TH2D* h2 = dH2(Form("%s/adcSpec_gem_%i_apv_%i",name, (int)i,
                (int)apv_id), Form("%s GEM Apv %i Spectrum",name,
                (int)apv_id), NUM_MAX_CHAN_PER_APV,SPEC_RANGE_MIN,
                NUM_MAX_CHAN_PER_APV, NUM_MAX_CHAN_PER_APV,RANGE_MIN,
                RANGE_MAX);
            temp_apv_hist_vec.push_back(h2);
        }
        adc_spectra.push_back(temp_apv_hist_vec);
        gem_spectra_x.push_back(gem_h2_x);
        gem_spectra_y.push_back(gem_h2_y);

        gem_spectra_x_cmode.push_back(gem_h2_x_cmode);
        gem_spectra_y_cmode.push_back(gem_h2_y_cmode);

        gem_hitmaps.push_back(gem_hitmap);
        gem_clustermaps.push_back(gem_clustermap);
        gem_cluster_multiplicity.push_back(gem_cluster_mult);

        temp_apv_hist_vec.clear();
    }
    return 0;
}

// prep IO
Long_t gemControl::startup(){
    initInput();
    initOutput();
    for (int i = 0; i<GEM_BANK_NUM; i++){
      binBlobReader[i] = new binaryBlobReader();
      std::cout << "binBlobReader: " << binBlobReader[i] << std::endl;
    }
    gemo = new LumiGEM();
    makeBranch("LumiGEMhits", (TObject **)&gemo);

    return 0;
}

// read in calibration vals should be done for every event in the both process loop!
Long_t gemControl::readCommonMode(){
  for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
      gem_list[i]->clearGemData();
      for(std::vector<apv*>::size_type j=0; j < gem_list[i]->apv_list.size(); j++){
             gem_list[i]->apv_list[j]->cmode = cmode_data[i][j];
      }
  }
  return ok;
}

// read in calibration vals
Long_t gemControl::readPedestal(){

    std::ifstream pedestal_data("gemControl_pedestal.dat");
    if(!pedestal_data){
        std::cout << "no pedestal found... wrong file?\n";
        exit(-1);
    }
    std::string pedestal_line;
    uint32_t gem_id,apv_id, apv_ch;
    float pedestal_data_val;

    while(std::getline(pedestal_data,pedestal_line)){
        std::istringstream line_data(pedestal_line);
        line_data >> gem_id >> apv_id >> apv_ch >> pedestal_data_val;
        if(apv_ch < NUM_MAX_CHAN_PER_APV){
          gem_list[gem_id]->apv_list[apv_id]->loadPedestal(apv_ch, pedestal_data_val);
          std::cout << "gem: " << gem_id << " apv: " << apv_id << " ch: " << apv_ch
              << " pedestal read in: " << pedestal_data_val << std::endl;
        }
    }

    return ok;
}



Long_t gemControl::determineCommonMode(){

    // clear all APV data from last event
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        gem_list[i]->clearGemData();
        for(std::vector<apv*>::size_type j=0;
            j < gem_list[i]->apv_list.size(); j++){
            gem_list[i]->apv_list[j]->clearApvData();
        }
    }

    // fill APVs with raw data from input tree
    size_t apv_data_read = readRawEventData();
    if(apv_data_read == 0){
      //  std::cout << "no APV data for event: " << event_number
      //      << std::endl;
        event_number++;
        return 0;
    }

    // clear cmode data for every event!
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        for(std::vector<apv*>::size_type j=0; j < gem_list[i]->apv_list.size(); j++){
                cmode_data[i][j] = 0;
                cmode_data_ct[i][j]= 0;
        }
    }

    // fills apv specific histos
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        gem_list[i]->setGemAxisData();
        for(std::vector<apv*>::size_type j=0; j < gem_list[i]->apv_list.size(); j++){

            // point and counter
            uint32_t* data_ptr =
                gem_list[i]->apv_list[j]->apv_channel_data;
            uint32_t first_valid_chan = gem_list[i]->\
                apv_list[j]->num_chan_to_skip;

            for(uint32_t k=first_valid_chan;k<NUM_MAX_CHAN_PER_APV;k++){
                cmode_data[i][j]+=*(data_ptr + inverse_mapping[k]);
                cmode_data_ct[i][j]++;
            }
            cmode_data[i][j] = cmode_data[i][j]/cmode_data_ct[i][j];
        }
    }


    event_number++;
    return 0;
}

Long_t gemControl::determinePedestal(){
    // clear all APV data from last event
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        gem_list[i]->clearGemData();
        for(std::vector<apv*>::size_type j=0;
            j < gem_list[i]->apv_list.size(); j++){
            gem_list[i]->apv_list[j]->clearApvData();
        }
    }

    // fill APVs with raw data from input tree
    size_t apv_data_read = readRawEventData();
    if(apv_data_read == 0){
      //  std::cout << "no APV data for event: " << event_number
      //      << std::endl;
        event_number++;
        return 0;
    }

    // fills apv specific histos
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        gem_list[i]->setGemAxisData();
        for(std::vector<apv*>::size_type j=0;
            j < gem_list[i]->apv_list.size(); j++){

            // point and counter
            uint32_t* data_ptr =
                gem_list[i]->apv_list[j]->apv_channel_data;
            uint32_t first_valid_chan = gem_list[i]->\
                apv_list[j]->num_chan_to_skip;

            for(uint32_t k=first_valid_chan;k<NUM_MAX_CHAN_PER_APV;k++){
                pedestal_data[i][j][k]+=*(data_ptr + k) - cmode_data[i][j];
                pedestal_data_ct[i][j][k]++;
            }
        }
    }

    return 0;
}




Long_t gemControl::process_avg(){

    // TODO clear output tree data holder
    // TODO visco stuff, not really interested in touching this
    gemo->hits.clear();//clear output tree
    // clear all APV data from last event
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        gem_list[i]->clearGemData();
        for(std::vector<apv*>::size_type j=0;
            j < gem_list[i]->apv_list.size(); j++){
            gem_list[i]->apv_list[j]->clearApvData();
        }
    }

    // fill APVs with raw data from input tree
    size_t apv_data_read = readRawEventData();

    if(apv_data_read == 0){
      //  std::cout << "no APV data for event: " << event_number
      //      << std::endl;
        event_number++;
        return 0;
    }

    // fills apv specific histos

    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        gem_list[i]->setGemAxisData();
        for(std::vector<apv*>::size_type j=0;
            j < gem_list[i]->apv_list.size(); j++){
            // point and counter
            uint32_t* data_ptr =
                gem_list[i]->apv_list[j]->apv_channel_data;
            uint32_t first_valid_chan = gem_list[i]->\
                apv_list[j]->num_chan_to_skip;
            uint32_t cmode = gem_list[i]->apv_list[j]->cmode;
            for(uint32_t k=first_valid_chan;k<NUM_MAX_CHAN_PER_APV;k++){
              adc_spectra[i][j]->Fill(k, *(data_ptr+inverse_mapping[k]));
              adc_spectra[i][j]->Fill(k, (int)(*(data_ptr+inverse_mapping[k]) - cmode));
            }
        }
    }

    // now I can do gems separately
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        for(uint32_t j=0; j<NUM_CHAN_PER_AXIS; j++){

            uint32_t* x_ptr = gem_list[i]->getDataPointerForChannel(j,
                Axis::X);
            uint32_t x_cmode=gem_list[i]->getCmodeForChannel(j,Axis::X);
            uint32_t x_pedestal=gem_list[i]->getPedestalForChannel(j,Axis::X);

            uint32_t* y_ptr = gem_list[i]->getDataPointerForChannel(j,
                Axis::Y);
            uint32_t y_cmode=gem_list[i]->getCmodeForChannel(j,Axis::Y);
            uint32_t y_pedestal=gem_list[i]->getPedestalForChannel(j,Axis::Y);


            gem_spectra_x[i]->Fill(j,*x_ptr);
            gem_spectra_y[i]->Fill(j,*y_ptr);

            gem_spectra_x_cmode[i]->Fill(j,(int)(*x_ptr-x_cmode-x_pedestal));
            gem_spectra_y_cmode[i]->Fill(j,(int)(*y_ptr-y_cmode-y_pedestal));

        }
    }

    // now find clusters
    // generate hitmaps and we're done
    debug(1, "Event NUmber = %i \n", event_number);
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){
        gem_list[i]->findClusters(event_number);

        if(gem_list[i]->gem_clusters.size() == 0){
          debug(1, "Event =%i Cluster Size(%i) = %i \n",\
                event_number, i, gem_list[i]->gem_clusters.size());
          //return 0;
        }

        for(size_t j=0; j<gem_list[i]->gem_clusters.size(); j++){
            gem_cluster cluster = gem_list[i]->gem_clusters[j];
            debug(1,"GEM_ID = %i => x = %zu y = %zu\n",\
              cluster.gem_id, cluster.x_chan_num, cluster.y_chan_num);
        }
        gem_cluster_multiplicity[i]->Fill(gem_list[i]->gem_clusters.size());

        //output tree will contain 3 highest quality cluster or less if
        //gem_clusters.size() < 3;
        size_t max_cluster_num = (gem_list[i]->gem_clusters.size() < 3) \
              ? gem_list[i]->gem_clusters.size() : 3;

        for(size_t j = 0; j < max_cluster_num; j++){
          gem_cluster biggest_cluster = gem_list[i]->gem_clusters[j];
          gem_hitmaps[i]->Fill(biggest_cluster.x_chan_num, biggest_cluster.y_chan_num);
          gem_clustermaps[i]->Fill(-1.0*biggest_cluster.x_chan_num*0.4+50,-1.0*biggest_cluster.y_chan_num*0.4+50);

          LumiGEMhit ahit;
            ahit.GEMid = biggest_cluster.gem_id;
            ahit.xl = biggest_cluster.x_chan_num;
            ahit.yl = biggest_cluster.y_chan_num;
            ahit.xlerr = biggest_cluster.x_chan_std;
            ahit.ylerr = biggest_cluster.y_chan_std;
            ahit.quality = biggest_cluster.adc_xy_avg_weight;
            ahit.charge = -1;//to be determined
            ahit.ampl = -1;
            ahit.sigma =-1;
          gemo->hits.push_back(ahit);
        }
    }

    // TODO writeout the gem clusters
    event_number++;
    return 0;
}

// debugging things can go here as well mostly just writes out common mode
// values
Long_t gemControl::finalize(){
    //dumpAllHistos();
    dumpAllGemInfo();
    //dumpCmodeInfo();
    return 0;
}

Long_t gemControl::cmdline(char *cmd){
    return 0; // 0 = all ok
}

// wrapper for gem addAPV to get hooked in init script
Long_t gemControl::associateApvToGem(int trash, uint32_t gem_id,
    uint32_t apv_id, Axis apv_axis, uint32_t chan_to_skip){

    // check this is a valid gem_id
    assert(gem_id < NUM_GEMS*GEM_BANK_NUM);
    assert(gem_list[gem_id]!=NULL);

    // add apv to the gem
    int ret = gem_list[gem_id]->addApv(apv_axis, apv_id, chan_to_skip);
    if(ret){
        std::cout << "failed to add apv: " << apv_id << "of axis: "
            << apv_axis << "to gem: " << gem_id << std::endl;
        exit(-1);
    }
    return 0;
}

// IO
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~ */
//TODO figure out how to error check this root IO garbage
Long_t gemControl::initInput(){
  printf("InitInput!!!!!\n");
  for(int i=0; i< GEM_BANK_NUM; i++){ //Ievgen: loop over banks
    getBranchObject(GEM_bank[i], (TObject**)&binInputDataHook[i]);
    //assert(binInputDataHook[i]!=NULL);
    if(binInputDataHook[i]!=NULL){
      std::cout << GEM_bank[i]<<": "<<"binInputDataHook: " << binInputDataHook[i] << std::endl;
      runInfo[i] = (MRTRunInfo*)getFileObject("RunInfo");
      assert(runInfo[i] != NULL);
      std::cout << GEM_bank[i]<<": " <<"runInfo: " << runInfo[i] << std::endl;
    } else{
      printf("ERROR: in gemControl::initInput(). NO data in %s\n",GEM_bank[i]);
    }
  }
  printf("End InitInput!!!\n");
    return 0;
}

// TODO setup output tree for hits
Long_t gemControl::initOutput(){
    return 0;
}

// if dumping the raw data gets read out to a file, but the file gets
// appended so make sure to delete it between runs for a fresh set of data

//Ievgen: This function should be fixed in the future to handle multiple GEMs in one bank;
size_t gemControl::readRawEventData(bool dump){

    size_t data_read_size = 0;
/*
    for(int i=0; i<GEM_BANK_NUM; i++){
        //assert(binBlobReader[i] != NULL);
        //assert(binInputDataHook[i] != NULL);
        if(binInputDataHook[i]!=NULL && binBlobReader[i] != NULL){
          data_read_size += binBlobReader[i]->readBank(binInputDataHook[i]);
        } else{
          printf("ERROR in readRawEventData!\n");
        }
    }

    if(!data_read_size){
        return data_read_size;
    }
*/
    for(std::vector<gem*>::size_type i=0; i<gem_list.size(); i++){

        if(binInputDataHook[i]!=NULL && binBlobReader[i] != NULL){

          data_read_size += binBlobReader[i]->readBank(binInputDataHook[i]);

          for(std::vector<apv*>::size_type j=0;
            j < gem_list[i]->apv_list.size(); j++){

            uint32_t apv_id = gem_list[i]->apv_list[j]->id;
            gem_list[i]->apv_list[j]->loadApvData(binBlobReader[i]-> //we need to figure out how to fix it in the future for more that one GEM in bank
                rawDataBlob.apv_data[apv_id]);

            if(dump){
                std::ofstream data_outfile;
                data_outfile.open("gemControl_AllData.dat",
                    std::fstream::out | std::fstream::app);
                data_outfile << "gem: " << i << " apv: "
                    << apv_id << std::endl;
                for(int k=0; k<128; k++){
                    data_outfile << "raw data: "
                    <<  *(binBlobReader[i]->rawDataBlob.apv_data[apv_id] + k)
                        << std::endl;
                }
            }
          }
        } //close if statement;
    }
//    printf("END of gemControl!\n");
    return data_read_size;
}


// debugging tools
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void gemControl::dumpAllGemInfo(){
    for(std::vector<gem*>::size_type i=0; i < NUM_GEMS*GEM_BANK_NUM; i++){
        if(gem_list[i]){
            gem_list[i]->dumpGemInfo();
        }
    }
}

// Idk how to get a TH name
void gemControl::dumpAllHistos(){
    std::cout << "dumping all histos:\n";
    std::cout << "~~~~~~~~~~~~~~~~~~~\n";

    std::cout << "gem_spectra x...\n";
    for(std::vector<TH1D*>::size_type i=0; i<gem_spectra_x.size(); i++){
        std::cout << "got histo: " << gem_spectra_x[i] << std::endl;
    }
    std::cout << "gem_spectra y...\n";
    for(std::vector<TH1D*>::size_type i=0; i<gem_spectra_y.size(); i++){
        std::cout << "got histo: " << gem_spectra_y[i] << std::endl;
    }

    std::cout << "gem_spectra x cmode ...\n";
    for(std::vector<TH1D*>::size_type i=0; i<gem_spectra_x_cmode.size(); i++){
        std::cout << "got histo: " << gem_spectra_x_cmode[i] << std::endl;
    }
    std::cout << "gem_spectra y cmode ...\n";
    for(std::vector<TH1D*>::size_type i=0; i<gem_spectra_y_cmode.size(); i++){
        std::cout << "got histo: " << gem_spectra_y_cmode[i] << std::endl;
    }

    std::cout << "apv_spectra...\n";
    //TODO get the typing right on this, idk size type for nested vectors
    // same with getting size of one slice of the array
    for(uint32_t i = 0; i < NUM_GEMS*GEM_BANK_NUM; i++){
        for(uint32_t j=0; j < NUM_APV_PER_GEM; j++){
            std::cout << "gem: " << i << " apv " << j << std::endl;
            std::cout << "got histo: " << adc_spectra[i][j] << std::endl;
        }
    }
}



void gemControl::dumpPedestalInfo(){
    std::cout << "dumping pedestal info:\n";
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::ofstream pedestal_outfile("gemControl_pedestal.dat");
    printf("Filling gemControl_pedestal.dat file! \n");

    for(uint32_t i=0; i < NUM_GEMS*GEM_BANK_NUM; i++){
      for(uint32_t j=0; j<NUM_APV_PER_GEM; j++){
        for(uint32_t k=0; k<NUM_MAX_CHAN_PER_APV; k++){
            double pedestal = pedestal_data[i][j][k]/pedestal_data_ct[i][j][k];
            if(pedestal != 1)
              //we are subtracting cmode at the end because right now it is not calculated on event to event basis
              //ideally cmode should be calculated for every event and subtracted from every event
            pedestal_outfile << i << "\t" << j  << "\t" << k << "\t" << pedestal
                << std::endl;
        }
      }
    }
    pedestal_outfile.close();
    printf("END of Filling gemControl_pedestal.dat file! \n");
}



extern "C"{
    Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p){
        return (Plugin *) new gemControl(in,out,inf_,outf_,p);
    }
}

ClassImp(gemControl);
