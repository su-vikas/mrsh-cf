#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/range/iterator_range.hpp>

#include "../lib/cuckoofilter/src/cuckoofilter.h"
#include "../header/config.h"
#include "../header/fnv.h"
#include "timing.h"
#include "../lib/cuckoofilter/src/hashutil.h"

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} // namespace 
using cuckoofilter::CuckooFilter;

class mrsh_CF{
    // for all the mrsh related functions
    size_t total_items; 
    std::uint64_t insert_errors = 0;
    std::uint64_t total_chunks = 0;

    public:
        CuckooFilter<char[], 32> filter;
         
        inline std::uint64_t getTotalChunks() {
            return this->total_chunks;
        }
        inline std::uint64_t getErrorChunks() {
            return this->insert_errors;
        }

        mrsh_CF(const size_t total_items);
        void hashBufferToCuckooFilter(unsigned char *byte_buffer, std::uint64_t bytes_read );
        void addPathToCuckooFilter(std::string file_path);
        unsigned char* getFileBuffer(std::string file_path);
        std::uint32_t roll_hashx(std::uint8_t c, std::uint8_t window[], std::uint32_t rhData[]);
        void addHashToCuckooFilter(std::uint64_t hashValue);
        void cuckooFilterContainsHash(std::uint64_t hashValue);
        void bufferInCuckooFilter(unsigned char *byte_buffer, uint64_t bytes_read, std::string);
        void pathInCuckooFilter(std::string path);
        void getCuckooFilterInfo();
        void copyBlockToBuffer(unsigned char* byte_buffer,char *buff, std::uint32_t start, std::uint32_t stop);
        std::uint64_t generateIndexTagHash(const char* buff,std::uint32_t size);
        void printBuff(char* buff, std::uint32_t size);
        void readFingerprint(std::string file);
        void writeFingerprint(std::string file);
        void compareFingerprints(const mrsh_CF obj);

};

//constructor
mrsh_CF::mrsh_CF(const size_t total_items): filter(total_items){
    this->total_items = total_items;
}

void
mrsh_CF::compareFingerprints( mrsh_CF obj){
    filter.CompareCF(obj.filter);
}
void
mrsh_CF::readFingerprint(std::string file){
    filter.readFingerprint(file);
}

void 
mrsh_CF::writeFingerprint(std::string file){
    filter.writeFingerprint(file);
}

// calculate rolling hash
std::uint32_t 
mrsh_CF::roll_hashx(std::uint8_t c, std::uint8_t window[], std::uint32_t rhData[]){
    rhData[2] -= rhData[1];
    rhData[2] += (ROLLING_WINDOW * c);

    rhData[1] += c;
    rhData[1] -= window[rhData[0] % ROLLING_WINDOW];

    window[rhData[0] % ROLLING_WINDOW] = c;
    rhData[0]++;

    /* The original spamsum AND'ed this value with 0xFFFFFFFF which
       in theory should have no effect. This AND has been removed
       for performance (jk) */
    rhData[3] = (rhData[3] << 5); //& 0xFFFFFFFF;
    rhData[3] ^= c;

    return rhData[1] + rhData[2] + rhData[3];
}

unsigned char* 
mrsh_CF::getFileBuffer(std::string file_path){
    //std::uint64_t bytes_read;       // stores the number of characters read from input file
    std::uint64_t file_size;        // store file size
    char *byte_buffer = NULL;   // TODO using signed char, check it doesnt break anything
    file_size = boost::filesystem::file_size(file_path);

    try{
        byte_buffer = new char[file_size+1];
    }
    catch (std::bad_alloc& ba){
        std::cerr << "allocation error for file size" << ba.what() << "\n";
    }

    std::ifstream is (file_path, std::ifstream::binary);
    if(is){
        is.read(byte_buffer, file_size);
        //bytes_read = is.gcount();
    }
    else{
        std::cout<< "Error in reading file \n";
    }
    return (unsigned char*)byte_buffer;
    //hashBufferToCuckooFilter(byte_buffer, bytes_read);
}

//TODO adapt to read any buffer, not just file
void 
mrsh_CF::hashBufferToCuckooFilter(unsigned char* byte_buffer, std::uint64_t bytes_read ){
    //TODO get file handle and size
    std::uint32_t i;                
    std::uint32_t last_block_index =0;
    std::uint64_t rValue;
    //std::uint64_t hashValue = 0;
    std::uint64_t chunkCount = 0;

    /* we need this arrays for extended rollhash function */
    std::uint8_t window[ROLLING_WINDOW] = {0};
    std::uint32_t rhData[4]             = {0};


    std::uint64_t hv;       // hash value
    std::uint32_t block_size; // size of the buffer block currently processed. 

    //TODO optimize the last block situation. code is being repeated unnecessariliy. 
    for (i=0; i<bytes_read; i++){
        rValue = roll_hashx(byte_buffer[i], window, rhData);
        if (rValue %BLOCK_SIZE == BLOCK_SIZE -1) {
            total_chunks++;
            block_size = i - last_block_index;

            //TODO untested
            //std::cout<< "block_size" << block_size << std::endl;
            char *buff = new char[block_size+1];//  = new char[block_size]; //+1 as last byte of the block is included. 
            copyBlockToBuffer(byte_buffer, buff, last_block_index, i);
            hv = generateIndexTagHash(buff, block_size);
            //std::cout<< "TagHash: " << hv <<std::endl;
            if(filter.Add(hv) != cuckoofilter::Ok){
                ++insert_errors;
                //std::cout << "\n ------------------------" << std::endl;
                //std::cout<<"Error in adding hash to cuckoo filter \n";
                //printBuff(buff, block_size);
            }
            //set indicies new
            last_block_index = i+1;
            if(i+SKIPPED_BYTES < bytes_read)
                i += SKIPPED_BYTES;

            delete[] buff;
        }
        hv = 0;
        block_size = 0;
    }
    
    //for the last block :-/
    //std::cout<<"last block index:"<<last_block_index<<"\t index:"<<i<< "\n"; 
    block_size = i - last_block_index;
    char *buff = new char[block_size+1]; 

    copyBlockToBuffer(byte_buffer, buff,last_block_index, i);
    hv = generateIndexTagHash(buff, block_size);
    //std::cout<< "TagHash: " << hv <<std::endl;
    if(filter.Add(hv) != cuckoofilter::Ok){
        ++insert_errors;
        //std::cout << " \n------------------------" << std::endl;
        //std::cout<<"\n Error in adding hash to cuckoo filter \n";
        //printBuff(buff, block_size);

    }

    //std::cout << "Total Chunks:"<<total_chunks << "\n";
    //std::cout << "Insert Errors:"<<insert_errors << "\n";
    delete[] buff;
    delete[] byte_buffer;
}

void
mrsh_CF::printBuff(char* buff, std::uint32_t size){
    std::uint32_t i;
    for(i=0;i<size;i++){
        std::cout<<buff[i];
    }
    std::cout<<std::endl;
}

void
mrsh_CF::copyBlockToBuffer(unsigned char* byte_buffer, char *buff, std::uint32_t start, std::uint32_t stop){
    // make the char[] of same size as the block size, so that sizeof() will give correct value. 
    // apparently not working with templates as char[] is unexpected
    std::uint32_t j;
    for(j=start;j<=stop;j++){
        buff[j-start] = byte_buffer[j];
    }
}

// generate the tag and index for a char buff, using SHA1 hashing
std::uint64_t 
mrsh_CF::generateIndexTagHash(const char* buff, std::uint32_t size) {
    //std::cout<<buff<<std::endl;
    std::string hashed_key = cuckoofilter::HashUtil::SHA1Hash((const char*) buff, size);
    //TODO check the pointer stuff over here with the defualt implementation
    std::uint64_t hv;
    hv = *((std::uint64_t*) hashed_key.c_str());
    //std::cout<< "Hash Value: " << hv <<std::endl;
    return hv;
}
// to add a file to CF, when path given
void 
mrsh_CF::addPathToCuckooFilter(std::string file_path){
    if (boost::filesystem::is_directory(file_path)){
        for (auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(file_path), {})) {
            // if a file, add files to the CF
            std::string file_str = entry.path().string();
            //std::cout<<"File: " << file_str << std::endl;
            unsigned char *byte_buffer;
            uint64_t file_size;

            byte_buffer = getFileBuffer(file_str);
            file_size = boost::filesystem::file_size(file_str);
            hashBufferToCuckooFilter(byte_buffer, file_size);
        }
    }
    else{
        unsigned char *byte_buffer;
        uint64_t file_size;

        byte_buffer = getFileBuffer(file_path);
        file_size = boost::filesystem::file_size(file_path);
        hashBufferToCuckooFilter(byte_buffer, file_size);
    }
    
}

void
mrsh_CF::pathInCuckooFilter(std::string path){
    if (boost::filesystem::is_directory(path)){
        for (auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path), {})) {
            std::string file_str = entry.path().string();
            unsigned char *byte_buffer;
            uint64_t file_size;

            byte_buffer = getFileBuffer(file_str);
            file_size = boost::filesystem::file_size(file_str);
            bufferInCuckooFilter(byte_buffer, file_size, file_str);
        }
    }
    else{
        unsigned char *byte_buffer;
        uint64_t file_size;

        byte_buffer = getFileBuffer(path);
        file_size = boost::filesystem::file_size(path);
        bufferInCuckooFilter(byte_buffer, file_size, path);
    }
}

std::uint64_t 
getDirSize(std::string path){
    std::uint64_t dir_size = 0;
    if (boost::filesystem::is_directory(path)){
        for (auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path), {})) {
            std::string file_str = entry.path().string();
            std::uint64_t file_size;

            file_size = boost::filesystem::file_size(file_str);
            dir_size += file_size;
        }
    }
    else{
        dir_size = boost::filesystem::file_size(path);
    }
    return dir_size;
}

/* CUCKOO FILTER RELATED CODE */
/*
void
mrsh_CF::addHashToCuckooFilter(std::uint64_t hashValue){
    if(filter.Add(hashValue) != cuckoofilter::Ok){
        std::cout<<"Error in adding hash to cuckoo filter \n";
    }
}
*/

/* check if cuckoo filter contains a hash */
void 
mrsh_CF::bufferInCuckooFilter(unsigned char *byte_buffer, uint64_t bytes_read, std::string filename){
    std::uint32_t chunkCount = 0; //how many hash value of chuncks present. 
    std::uint32_t chunkDetected = 0; //how many hash value of chuncks present. 

    std::uint32_t i;                
    std::uint32_t last_block_index =0;
    std::uint64_t rValue, hashValue = 0;

    /* we need this arrays for extended rollhash function */
    std::uint8_t window[ROLLING_WINDOW] = {0};
    std::uint32_t rhData[4]             = {0};

    std::uint32_t block_size;
    std::uint64_t hv;

    for (i=0; i<bytes_read; i++){
        rValue = roll_hashx(byte_buffer[i], window, rhData);
        if (rValue %BLOCK_SIZE ==  BLOCK_SIZE -1) {
            std::uint32_t j;

            block_size = i - last_block_index;
            char *buff = new char[block_size+1];
            copyBlockToBuffer(byte_buffer, buff, last_block_index, i);

            ++chunkCount;
            hv = generateIndexTagHash(buff, block_size);
            if(filter.Contain(hv) == cuckoofilter::Ok){
                chunkDetected++;
            }
            else{
                //std::cout<<hv;
                //std::cout<<chunkCount<<"\t" <<buff << "\n";
            }
            
            //set indicies new
            last_block_index = i+1;
            if(i+SKIPPED_BYTES < bytes_read)
                i += SKIPPED_BYTES;
            delete[] buff;
        }
        hv = 0;
        block_size =0;
    }
    //FOR LAST BLOCK

    block_size = i - last_block_index;
    char *buff = new char[block_size+1]; 

    copyBlockToBuffer(byte_buffer, buff,last_block_index, i);
    ++chunkCount;
    hv = generateIndexTagHash(buff, block_size);
    if(filter.Contain(hv) == cuckoofilter::Ok){
        chunkDetected++;
    }
    else{
        //std::cout<<hv;
        //std::cout<<chunkCount<<"\t" <<buff << "\n";
    }


    //if(chunkCount+1 != chunkDetected){
//std::cout<< filename<<  "\t Total Chunks: "<< chunkCount <<  "\t Chunks Detected: " << chunkDetected<<"\n";
    std::cout<< filename<<  ","<< chunkCount <<  "," << chunkDetected<<"\n";
    //}

    //std::cout<< "Total Chunks: "<< chunkCount+1 <<  "\t Chunks Detected: " << chunkDetected<<"\n";
    delete[] buff;
    delete[] byte_buffer;
}

void 
mrsh_CF::getCuckooFilterInfo(){
    std::cout<< "Load Factor is: "<< filter.Info() << "\n";
}

int main(int argc, char **argv){
    try{
        namespace po = boost::program_options;
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "Print help message")
            ("file,f", po::value<std::string>()->required(),"Set of files/dir against which a file can be compared")
            ("signature,s", po::value<std::string>(),"Give first signature file")
            ("signature2,r", po::value<std::string>(),"Give second signature file")
            ("generate,g", po::value<std::string>(), "Generate signature and write to a file for file given with -f")
            ("info,i", "Show cuckoo filter related info")
            ("compare,c", po::value<std::string>()->required(),"file/dir to be compared against file given with -f");

        po::variables_map vm;
        try{
            po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
        }

        catch(boost::program_options::required_option& e) 
        { 
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
            return ERROR_IN_COMMAND_LINE; 
        } 
        catch(boost::program_options::error& e) 
        { 
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
            return ERROR_IN_COMMAND_LINE; 
        }

        if (vm.count("help")) {
            std::cout << "MRSH with Cuckoo Filter\n";
            std::cout << desc;
            return 0;
        }
        //std::cout<< "Files: " << vm["compare"].as<std::string>()<< std::endl;
        std::uint64_t dir_size;

        if (vm.count("file")){
            dir_size = getDirSize(vm["file"].as<std::string>());
        }

        mrsh_CF obj(dir_size/BLOCK_SIZE/2);
        mrsh_CF obj2(dir_size/BLOCK_SIZE/2);
        //mrsh_CF obj(3000);
        //mrsh_CF obj2(3000);
        //std::cout << (dir_size/BLOCK_SIZE/2) << std::endl;
        if (vm.count("file")){
            obj.addPathToCuckooFilter(vm["file"].as<std::string>());
            //obj2.addPathToCuckooFilter(vm["compare"].as<std::string>());
            //obj.compareFingerprints(obj2);
            /*std::cout<<vm["file"].as<std::string>()<<",";*/
        }

        if(vm.count("signature") && vm.count("signature2")){
            obj.readFingerprint(vm["signature"].as<std::string>());
            obj2.readFingerprint(vm["signature2"].as<std::string>());
            obj.compareFingerprints(obj2);
        }

        if (vm.count("generate")){
            obj.writeFingerprint(vm["generate"].as<std::string>());
        }

        if (vm.count("info")){
            obj.getCuckooFilterInfo();
            std::cout<<"total chunks:" << obj.getTotalChunks() << "\t errors: "<<obj.getErrorChunks() <<std::endl; 
        }

        if(vm.count("compare")){
            obj.pathInCuckooFilter(vm["compare"].as<std::string>());
        }
    }
    catch(std::exception& e)
    {
        std::cerr<<"[!] Error occurred "<<e.what() << std::endl;
    }

}
