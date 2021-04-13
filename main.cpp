//
//  main.cpp
//  GNQTS RSI
//
//  Created by 唐健恆 on 2021/4/9.
//  Copyright © 2021 Alvin. All rights reserved.
//

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cfloat>
#include "date.h"

using namespace std;
using namespace filesystem;

#define EXPNUMBER 50
#define ITERNUMBER 10000
#define PARTICLENUMBER 10
#define FUNDS 10000000
#define DELTA 0.004
#define QTSTYPE 2 //QTS 0, GQTS 1, GNQTS 2
//#define TRENDLINETYPE 0 //linear 0, quadratic 1
#define STARTDATE "2010-01-04"
#define RSI_RANGE 256
#define RSI_BIT_SIZE 8
#define BUY_BIT_SIZE 8
#define SELL_BIT_SIZE 8
#define MODE "exhaustive" //train, exhaustive

string FILE_DIR = "0412_sliding_window_data_temp";
string COMPANY_PRICE_DIR = "select_stock_price";
string RSI_DIR = "select_RSI_list";
int BIT_SIZE = RSI_BIT_SIZE + BUY_BIT_SIZE + SELL_BIT_SIZE;


class CompanyData{
public:
    string *date_list = NULL;
    double *price_list = NULL;
    double **RSI_list = NULL;
    int day_number;
    void init(int);
    void init(CompanyData&);
    void print();
    
    ~CompanyData();
};

CompanyData::~CompanyData(){
    delete[] this -> date_list;
    delete[] this -> price_list;
    for(int j = 0; j < RSI_RANGE; j++){
        delete[] this -> RSI_list[j];
    }
    delete[] this -> RSI_list;
}

void CompanyData::init(int day_number){
    this -> day_number = day_number;
    this -> date_list = new string[day_number];
    this -> price_list = new double[day_number];
    this -> RSI_list = new double*[RSI_RANGE];
    for(int j = 0; j < RSI_RANGE; j++){
        RSI_list[j] = new double[day_number];
    }
}

void CompanyData::init(CompanyData &c){
    this -> init(c.day_number);
    this -> day_number = c.day_number;
    for(int j = 0; j < this -> day_number; j++){
        this -> date_list[j] = c.date_list[j];
        this -> price_list[j] = c.price_list[j];
    }
    
    for(int j = 0; j < RSI_RANGE; j++){
        for(int k = 0; k < this -> day_number; k++){
            this -> RSI_list[j][k] = c.RSI_list[j][k];
        }
    }
}

void CompanyData::print(){
    cout << "day_number: " << day_number << endl;
    for(int j = 0; j < day_number; j++){
        cout << date_list[j] << ", " << price_list[j] << ", " << RSI_list[13][j] << endl;
    }
}

class RSIParticle{
public:
    int gen = 0;
    int exp = 0;
    int answer_counter = 0;
    int day_number = 0;
    int trade_times = 0;
    int RSI_number = 0;
    int investment_number = 0;
    int bit_size = 0;
    int* investment_list = NULL;
    int* data = NULL;
    int* trade_list = NULL;
    int* trade_record = NULL;
    double funds = 0;
    double return_rate = 0;
    double upper_bound = 0;
    double lower_bound = 0;
    double remain_money = 0;
    double* remain_money_list = NULL;
    double* total_money = NULL;
    
    CompanyData companyData;
    
    
    
    RSIParticle();
    ~RSIParticle();
    RSIParticle(int, int, double, CompanyData&);
    void init();
    void init(int, int, double, CompanyData&);
    void bitToDec();
    void print();
    void copyP(RSIParticle&);
    double getRSI(int);
    double getReturnRate();
};


RSIParticle::RSIParticle(){
    
}

RSIParticle::RSIParticle(int day_number, int bit_size, double funds, CompanyData &companyData){
    
    this -> day_number = day_number;
    this -> funds = funds;
    this -> remain_money = funds;
    this -> bit_size = bit_size;
    if(this -> data != NULL){
        delete[] this -> data;
        delete[] this -> trade_list;
        delete[] this -> investment_list;
        delete[] this -> total_money;
        delete[] this -> remain_money_list;
    }
    this -> data = new int[bit_size];
    this -> trade_list = new int[day_number];
    this -> investment_list = new int[day_number];
    this -> total_money = new double[day_number];
    this -> remain_money_list = new double[day_number];
    this -> companyData.init(companyData);
}

RSIParticle::~RSIParticle(){
    if(this -> data != NULL){
        delete[] this -> data;
        delete[] this -> trade_list;
        delete[] this -> investment_list;
        delete[] this -> total_money;
        delete[] this -> remain_money_list;
    }
    this -> data = NULL;
    this -> trade_list = NULL;
    this -> investment_list = NULL;
    this -> total_money = NULL;
    this -> remain_money_list = NULL;
    
}

void RSIParticle::init(){
    this -> gen = 0;;
    this -> exp = 0;;
    this -> answer_counter = 0;
    this -> trade_times = 0;
    this -> RSI_number = 0;
    this -> investment_number = 0;
    this -> return_rate = 0;
    this -> upper_bound = 0;
    this -> lower_bound = 0;
    this -> remain_money =  this -> funds;
    if(this -> data != NULL){
        delete[] this -> data;
        delete[] this -> trade_list;
        delete[] this -> investment_list;
        delete[] this -> total_money;
        delete[] this -> remain_money_list;
    }
    this -> data = new int[this -> bit_size];
    this -> trade_list = new int[this -> day_number];
    this -> investment_list = new int[this -> day_number];
    this -> total_money = new double[this -> day_number];
    this -> remain_money_list = new double[this -> day_number];
    
}

void RSIParticle::init(int day_number, int bit_size, double funds, CompanyData &companyData){
    this -> gen = 0;;
    this -> exp = 0;;
    this -> answer_counter = 0;
    this -> trade_times = 0;
    this -> RSI_number = 0;
    this -> investment_number = 0;
    this -> return_rate = 0;
    this -> upper_bound = 0;
    this -> lower_bound = 0;
    this -> remain_money = funds;
    this -> day_number = day_number;
    this -> funds = funds;
    this -> companyData.init(companyData);
    this -> bit_size = bit_size;
    
    if(this -> data != NULL){
        delete[] this -> data;
        delete[] this -> trade_list;
        delete[] this -> investment_list;
        delete[] this -> total_money;
        delete[] this -> remain_money_list;
    }
    this -> data = new int[bit_size];
    this -> trade_list = new int[day_number];
    this -> investment_list = new int[day_number];
    this -> total_money = new double[day_number];
    this -> remain_money_list = new double[day_number];
}

void RSIParticle::bitToDec(){
    int shift = 0;
    int sum = 0;
    for(int j = 0; j < RSI_BIT_SIZE; j++){
        int temp = j + shift;
        if(this -> data[temp] == 1){
            sum += pow(2, RSI_BIT_SIZE - j - 1);
        }
    }
    this -> RSI_number = sum + 1;
    sum = 0;
    shift += RSI_BIT_SIZE;
    for(int j = 0; j < BUY_BIT_SIZE; j++){
        int temp = j + shift;
        if(this -> data[temp] == 1){
            sum += pow(2, BUY_BIT_SIZE - j - 1);
        }
    }
    this -> lower_bound = sum;
    sum = 0;
    shift += BUY_BIT_SIZE;
    for(int j = 0; j < SELL_BIT_SIZE; j++){
        int temp = j + shift;
        if(this -> data[temp] == 1){
            sum += pow(2, SELL_BIT_SIZE - j - 1);
        }
    }
    this -> upper_bound = sum;
}

void RSIParticle::print(){
    cout << "exp: " << this -> exp << endl;
    cout << "gen: " << this -> gen << endl;
    cout << "data: " << endl;
    for(int j = 0; j < BIT_SIZE; j++){
        if(j > 0){
            cout << ", ";
        }
        cout << this -> data[j];
    }
    cout << endl;
    cout << "RSI_number: " << this -> RSI_number << endl;
    cout << "LowerBound: " << this -> lower_bound << endl;
    cout << "UpperBound: " << this -> upper_bound << endl;
    cout << "Trade times: " << this -> trade_times << endl;
    cout << "Return rate: " << this -> return_rate << endl;
    cout << endl;
}

double RSIParticle::getRSI(int day_number){
    return this -> companyData.RSI_list[this -> RSI_number - 1][day_number];
}

double RSIParticle::getReturnRate(){
    return (this -> total_money[this -> day_number - 1] - this -> funds) / this -> funds * 100;
}

void RSIParticle::copyP(RSIParticle& a) {
    this -> gen = a.gen;
    this -> exp = a.exp;
    this -> answer_counter = a. answer_counter;
    this -> day_number = a.day_number;
    this -> trade_times = a.trade_times;
    this -> RSI_number = a.RSI_number;
    this -> investment_number = a.investment_number;
    this -> funds = a.funds;
    this -> return_rate = a.return_rate;
    this -> upper_bound = a.upper_bound;
    this -> lower_bound = a.lower_bound;
    this -> remain_money = a.remain_money;
    for(int j = 0; j < BIT_SIZE; j++){
        this -> data[j] = a.data[j];
    }
    for(int j = 0; j < this -> day_number; j++){
        this -> investment_list[j] = a.investment_list[j];
        this -> trade_list[j] = a.trade_list[j];
        this -> remain_money_list[j] = a.remain_money_list[j];
        this -> total_money[j] = a.total_money[j];
    }
}




string* vectorToArray(vector<string> &data_vector){
    string *d = new string[data_vector.size()];
    for(int j = 0; j < data_vector.size(); j++){
        d[j] = data_vector[j];
    }
    return d;
}

string** vectorToArray(vector<vector<string>> &data_vector){
    string **d = new string*[data_vector.size()];
    for(int j = 0; j < data_vector.size(); j++){
        d[j] = vectorToArray(data_vector[j]);
    }
    return d;
}

string*** vectorToArray(vector<vector<vector<string>>> &data_vector){
    string ***d = new string**[data_vector.size()];
    for(int j = 0; j < data_vector.size(); j++){
        d[j] = vectorToArray(data_vector[j]);
    }
    return d;
}

vector<string> genCompanyList(string fileDir){
    cout << fileDir << endl;
    vector<string> temp;
    path p1(fileDir);
    for( directory_iterator it = directory_iterator(p1);
        it != directory_iterator(); ++ it )
    {
        path px = it->path();
        if(px.extension() == ".csv"){
            temp.push_back(px.stem().generic_string());
        }
    }
    sort(temp.begin(), temp.end());
    return temp;
}

bool readData(string filename, vector<vector<string>> &data_vector) {
    cout << filename << endl;
    ifstream inFile(filename, ios::in);
    string line;
    vector< vector<string> > temp;

    if (!inFile) {
        cout << "Open file failed!" << endl;
        exit(1);
    }
    while (getline(inFile, line)) {
        istringstream delime(line);
        string s;
        vector<string> line_data;
        while (getline(delime, s, ',')) {
            if (s != "\r") {
                s.erase(remove(s.begin(), s.end(), '\r'), s.end());
                line_data.push_back(s);
            }
        }
        temp.push_back(line_data);
    }
    inFile.close();
    data_vector = temp;
    return true;
}

bool readRSIData(string file_dir, vector<vector<vector<string>>> &RSI_data_vector) {
    vector<vector<string>> data_vector;
    for(int j = 0; j < RSI_RANGE; j++){
        string file_name = file_dir + "RSI(" + to_string(j+1) + ").csv";
        if(!readData(file_name, data_vector)){
            cout << "Open file failed!" << endl;
            exit(1);
        }
        RSI_data_vector.push_back(data_vector);
        data_vector.clear();
    }
    return true;
}

void createCompanyData(CompanyData &companyData, string **price_data, string ***RSI_data, int start_index, int end_index){
    for(int j = start_index; j <= end_index ; j++){
        companyData.date_list[j - start_index] = price_data[j][0];
        companyData.price_list[j - start_index] = atof(price_data[j][4].c_str());
        int temp = j - 2;
        for(int k = 0; k < RSI_RANGE; k++){
            companyData.RSI_list[k][j - start_index] = atof(RSI_data[k][temp - k][1].c_str());
        }
    }
}

void initial(double *b, int size) {
    for (int j = 0; j < size; j++) {
        b[j] = 0.5;
    }
}

void initRSIParticle(RSIParticle *p, int particle_number, int day_number, CompanyData &companyData){
    for(int j = 0; j < particle_number; j++){
        p[j].init(day_number, BIT_SIZE, FUNDS, companyData);
    }
}

void initRISParticle(RSIParticle *p){
    for(int j = 0; j < PARTICLENUMBER; j++){
        p[j].init();
    }
}

void genParticle(RSIParticle *p, int n, int i, double *beta_){
    for (int j = 0; j < PARTICLENUMBER; j++) {
        p[j].exp = n + 1;
        p[j].gen = i + 1;
        for (int k = 0; k < BIT_SIZE; k++) {
            double r = (double)rand() / (double)RAND_MAX;
            if (r > beta_[k]) {
                p[j].data[k] = 0;
            }
            else {
                p[j].data[k] = 1;
            }
        }
        p[j].bitToDec();
    }
}

void startTrade(RSIParticle *p, int day_number, int particle_number){
    
    for(int j = 0; j < particle_number; j++){
        for(int k = 0; k < day_number; k++){
            if(p[j].getRSI(k) <= p[j].lower_bound && p[j].investment_number == 0 && k != day_number - 1){
                p[j].investment_number = p[j].remain_money / p[j].companyData.price_list[k];
                p[j].remain_money -= p[j].investment_number * p[j].companyData.price_list[k];
                p[j].trade_list[k] = 1;
                p[j].trade_times += 1;
            }else if(p[j].getRSI(k) >= p[j].upper_bound && p[j].investment_number != 0){
                p[j].remain_money += p[j].investment_number * p[j].companyData.price_list[k];
                p[j].investment_number = 0;
                p[j].trade_list[k] = 1;
            }else{
                p[j].trade_list[k] = 0;
            }
            
            if(k == day_number - 1 && p[j].investment_number != 0){
                p[j].remain_money += p[j].investment_number * p[j].companyData.price_list[k];
                p[j].investment_number = 0;
                p[j].trade_list[k] = 1;
            }
            
            p[j].remain_money_list[k] = p[j].remain_money;
            p[j].investment_list[k] = p[j].investment_number;
            p[j].total_money[k] = p[j].investment_number * p[j].companyData.price_list[k] + p[j].remain_money;
        }
        p[j].return_rate = (p[j].remain_money - p[j].funds) / p[j].funds * 100;
    }
}

void recordGAnswer(RSIParticle *p, RSIParticle &gBest, RSIParticle &gWorst, RSIParticle &pBest, RSIParticle &pWorst) {
    pBest.copyP(p[0]);
    pWorst.copyP(p[0]);
    for (int j = 0; j < PARTICLENUMBER; j++) {
        if (pBest.return_rate < p[j].return_rate) {
            pBest.copyP(p[j]);
        }
        if (pWorst.return_rate > p[j].return_rate) {
            pWorst.copyP(p[j]);
        }
    }
    
    if (gBest.return_rate < pBest.return_rate) {
            gBest.copyP(pBest);
    }
    
    if (gWorst.return_rate > pWorst.return_rate) {
        gWorst.copyP(pWorst);
    }
}

void adjBeta(RSIParticle& best, RSIParticle& worst, double *beta_) {
    for (int j = 0; j < BIT_SIZE; j++) {
        if (QTSTYPE == 2) {
            if (best.data[j] > worst.data[j]) {
                if (beta_[j] < 0.5) {
                    beta_[j] = 1 - beta_[j];
                }
                beta_[j] += DELTA;
            }
            else if (best.data[j] < worst.data[j]) {
                if (beta_[j] > 0.5) {
                    beta_[j] = 1 - beta_[j];
                }
                beta_[j] -= DELTA;
            }
        }
        else if (QTSTYPE == 1) {
            if (best.data[j] > worst.data[j]) {
                beta_[j] += DELTA;
            }
            else if (best.data[j] < worst.data[j]) {
                beta_[j] -= DELTA;
            }
        }
        else {
            if (best.data[j] > worst.data[j]) {
                beta_[j] += DELTA;
            }
            else if (best.data[j] < worst.data[j]) {
                beta_[j] -= DELTA;
            }
        }
    }
}

void recordExpAnswer(RSIParticle& expBest, RSIParticle& gBest) {
    if (expBest.return_rate < gBest.return_rate) {
        expBest.copyP(gBest);
        expBest.answer_counter = 1;
    }
    else if (expBest.return_rate == gBest.return_rate) {
        expBest.answer_counter++;
    }
}

void genTradeRecord(RSIParticle& expBest, int day_number){
    bool flag = true;
    int counter = 0;
    expBest.trade_record = new int[expBest.trade_times * 2];
    for(int j = 0; j < day_number; j++){
        if(expBest.trade_list[j] == 1){
            expBest.trade_record[counter] = j;
            counter++;
        }
    }
}

void outputFile(RSIParticle &p, string file_name) {
    ofstream outfile;
    outfile.open(file_name, ios::out);
    outfile << fixed << setprecision(10);
    outfile << "Iteration," << ITERNUMBER << endl;
    outfile << "Element number," << PARTICLENUMBER << endl;
    outfile << "Delta," << DELTA << endl;
    outfile << "Exp times," << EXPNUMBER << endl;
    outfile << "QTS TYPE," << QTSTYPE << endl;
    outfile << endl;
    
    outfile << "Init funds," << p.funds << endl;
    outfile << "Final funds," << p.total_money[p.day_number - 1] << endl;
    outfile << "Real award," << p.total_money[p.day_number - 1] - p.funds << endl;
    outfile << endl;
    
    outfile << "RSI number," << p.RSI_number << endl;
    outfile << "Buy point," << p.lower_bound << endl;
    outfile << "Sell point," << p.upper_bound << endl;
    outfile << "Trade times," << p.trade_times << endl;
    outfile << "Return rate," << p.return_rate << endl;
    outfile << endl;
    
    outfile << "Best experiment," << p.exp << endl;
    outfile << "Best generation," << p.gen << endl;
    outfile << "Best answer times," << p.answer_counter << endl;
    outfile << endl;
    
    outfile << "Trade record,Date,Price,RSI,Stock number,Remain money,Total money," << endl;
    for(int j = 0; j < p.trade_times * 2; j += 2){
        outfile << "Buy," << p.companyData.date_list[p.trade_record[j]] << "," << p.companyData.price_list[p.trade_record[j]] << "," << p.companyData.RSI_list[p.RSI_number - 1][p.trade_record[j]] << "," << p.investment_list[p.trade_record[j]] << "," << p.remain_money_list[p.trade_record[j]] << "," << p.total_money[p.trade_record[j]] << "," << endl;
        outfile << "Sell," << p.companyData.date_list[p.trade_record[j+1]] << "," << p.companyData.price_list[p.trade_record[j+1]] << "," << p.companyData.RSI_list[p.RSI_number - 1][p.trade_record[j+1]] << "," << p.investment_list[p.trade_record[j+1]] << "," << p.remain_money_list[p.trade_record[j+1]] << "," << p.total_money[p.trade_record[j+1]] << "," << endl;
        outfile << endl;
    }
    outfile.close();
}


void recordCPUTime(double START, double END, string file_name){
    double total_time = (END - START) / CLOCKS_PER_SEC;
    ofstream outfile_time;
    outfile_time.open(file_name, ios::out);
    outfile_time << "total time: " << total_time << " sec" << endl;
    outfile_time.close();
}

void recordData(RSIParticle &expBest, ofstream &outfile_data){
    outfile_data << expBest.companyData.date_list[0] << " ~ ";
    outfile_data << expBest.companyData.date_list[expBest.day_number - 1] << ",";
    outfile_data << expBest.exp << "," << expBest.gen << ",";
    outfile_data << expBest.RSI_number << "," << expBest.lower_bound << "," << expBest.upper_bound << ",";
    outfile_data << expBest.trade_times << "," << expBest.return_rate << "," << endl;
}

void createDir(string file_dir, string company_name, string type){
    create_directory(file_dir);
    create_directory(file_dir + "/" + company_name);
    create_directory(file_dir + "/" + company_name + "/" + type);
    create_directory(file_dir + "/" + company_name + "/" + type + "/" + "train");
    create_directory(file_dir + "/" + company_name + "/" + type + "/" + "test");
    create_directory(file_dir + "/" + company_name + "/" + type + "/" + "exhaustive");
}

void preSet(string mode, Date& current_date, Date& finish_date, int SLIDETYPE, string& TYPE) {
    string STARTYEAR;
    string STARTMONTH;
    string ENDYEAR;
    string ENDMONTH;
    int slide_number;
    int train_range;
    switch (SLIDETYPE) {
        case 0:
            STARTYEAR = "2010";
            STARTMONTH = "1";
            ENDYEAR = "2020";
            ENDMONTH = "12";
            TYPE = "M2M";
            train_range = 1;
            slide_number = 1;
            break;

        case 1:
            STARTYEAR = "2010";
            STARTMONTH = "1";
            ENDYEAR = "2020";
            ENDMONTH = "12";
            TYPE = "Q2Q";
            train_range = 3;
            slide_number = 3;
            break;
            
        case 2:
            STARTYEAR = "2010";
            STARTMONTH = "1";
            ENDYEAR = "2020";
            ENDMONTH = "12";
            TYPE = "H2H";
            train_range = 6;
            slide_number = 6;
            break;
        
        case 3:
            STARTYEAR = "2010";
            STARTMONTH = "1";
            ENDYEAR = "2020";
            ENDMONTH = "12";
            TYPE = "Y2Y";
            train_range = 12;
            slide_number = 12;
            break;
            
        case 4:
        STARTYEAR = "2010";
        STARTMONTH = "1";
        ENDYEAR = "2020";
        ENDMONTH = "12";
        TYPE = "A2A";
        train_range = 132;
        slide_number = 132;
        break;
        
    }
//    switch (SLIDETYPE) {
//    case 0:
//        STARTYEAR = "2009";
//        STARTMONTH = "12";
//        ENDYEAR = "2019";
//        ENDMONTH = "11";
//        TYPE = "M2M";
//        train_range = 1;
//        slide_number = 1;
//        break;
//    case 1:
//        STARTYEAR = "2009";
//        STARTMONTH = "10";
//        ENDYEAR = "2019";
//        ENDMONTH = "9";
//        TYPE = "Q2M";
//        train_range = 3;
//        slide_number = 1;
//        break;
//    case 2:
//        STARTYEAR = "2009";
//        STARTMONTH = "10";
//        ENDYEAR = "2019";
//        ENDMONTH = "7";
//        TYPE = "Q2Q";
//        train_range = 3;
//        slide_number = 3;
//        break;
//    case 3:
//        STARTYEAR = "2009";
//        STARTMONTH = "7";
//        ENDYEAR = "2019";
//        ENDMONTH = "6";
//        TYPE = "H2M";
//        train_range = 6;
//        slide_number = 1;
//        break;
//    case 4:
//        STARTYEAR = "2009";
//        STARTMONTH = "7";
//        ENDYEAR = "2019";
//        ENDMONTH = "4";
//        TYPE = "H2Q";
//        train_range = 6;
//        slide_number = 3;
//        break;
//    case 5:
//        STARTYEAR = "2009";
//        STARTMONTH = "7";
//        ENDYEAR = "2019";
//        ENDMONTH = "1";
//        TYPE = "H2H";
//        train_range = 6;
//        slide_number = 6;
//        break;
//    case 6:
//        STARTYEAR = "2009";
//        STARTMONTH = "1";
//        ENDYEAR = "2018";
//        ENDMONTH = "12";
//        TYPE = "Y2M";
//        train_range = 12;
//        slide_number = 1;
//        break;
//    case 7:
//        STARTYEAR = "2009";
//        STARTMONTH = "1";
//        ENDYEAR = "2018";
//        ENDMONTH = "10";
//        TYPE = "Y2Q";
//        train_range = 12;
//        slide_number = 3;
//        break;
//    case 8:
//        STARTYEAR = "2009";
//        STARTMONTH = "1";
//        ENDYEAR = "2018";
//        ENDMONTH = "7";
//        TYPE = "Y2H";
//        train_range = 12;
//        slide_number = 6;
//        break;
//    case 9:
//        STARTYEAR = "2009";
//        STARTMONTH = "1";
//        ENDYEAR = "2018";
//        ENDMONTH = "1";
//        TYPE = "Y2Y";
//        train_range = 12;
//        slide_number = 12;
//        break;
//    case 10:
//        STARTYEAR = "2009";
//        STARTMONTH = "1";
//        ENDYEAR = "2018";
//        ENDMONTH = "12";
//        TYPE = "M#";
//        if(mode == "train"){
//            train_range = 1;
//        }else{
//            train_range = 12;
//        }
//        slide_number = 1;
//        break;
//    case 11:
//        STARTYEAR = "2009";
//        STARTMONTH = "1";
//        ENDYEAR = "2018";
//        ENDMONTH = "10";
//        TYPE = "Q#";
//        if(mode == "train"){
//            train_range = 3;
//        }else{
//            train_range = 12;
//        }
//        slide_number = 3;
//        break;
//    case 12:
//        STARTYEAR = "2009";
//        STARTMONTH = "1";
//        ENDYEAR = "2018";
//        ENDMONTH = "7";
//        TYPE = "H#";
//        if(mode == "train"){
//            train_range = 6;
//        }else{
//            train_range = 12;
//        }
//        slide_number = 6;
//        break;
//    }

    current_date.date.tm_year = atoi(STARTYEAR.c_str()) - 1900;
    current_date.date.tm_mon = atoi(STARTMONTH.c_str()) - 1;
    current_date.date.tm_mday = 1;
    current_date.slide_number = slide_number;
    current_date.train_range = train_range;
    
    finish_date.date.tm_year = atoi(ENDYEAR.c_str()) - 1900;
    finish_date.date.tm_mon = atoi(ENDMONTH.c_str()) - 1;
    finish_date.date.tm_mday = 1;
    finish_date.slide_number = slide_number;
    finish_date.train_range = train_range;
    
    if(mode == "test"){
        current_date.slide(train_range);
        finish_date.slide(train_range);
    }
}

void setWindow(string mode, string &start_date, string &end_date, int &start_index, int &end_index, Date current_date, Date finish_data, string* data_copy, string **data, int day_number, int &range_day_number){
    bool sw = false;//判斷是否找到開始日期
    int flag = 0;//判斷是否找到結束月份
    for(int j = 0; j < day_number; j++){
        string temp1 = current_date.getYear() + "-" + current_date.getMon();
        string temp2;
        if(mode == "train"){
            temp2 = current_date.getRangeEnd(current_date.train_range - 1).getYear() + "-" + current_date.getRangeEnd(current_date.train_range - 1).getMon();
        }else{
            temp2 = current_date.getRangeEnd(current_date.slide_number - 1).getYear() + "-" + current_date.getRangeEnd(current_date.slide_number - 1).getMon();
        }
        
        if(!sw && temp1 == data_copy[j]){
            start_date = data[j+1][0];
            start_index = j + 1;
            sw = true;
        }
        
        if(sw){
            if(flag == 1 && temp2 != data_copy[j]){
                end_date = data[j][0];
                end_index = j;
                flag = 0;
                break;
            }
            if(flag == 0 && temp2 == data_copy[j]){
                flag = 1;
            }
        }
    }
    
    if(flag == 1){
        end_date = data[day_number][0];
        end_index = day_number;
    }
    range_day_number = end_index - start_index + 1;
}

void copyData(string *data_copy, string **data, int day_number){
    for(int j = 0; j < day_number; j++){
        data_copy[j] = data[j+1][0];
        data_copy[j].resize(7);
    }
}

void releaseData(vector<vector<string>> &price_data_vector, vector<vector<vector<string>>> &RSI_data_vector, string **price_data, string ***RSI_data, string *data_copy){
    for(int j = 0; j < price_data_vector.size(); j++){
        delete[] price_data[j];
        price_data_vector[j].clear();
    }
    delete[] price_data;
    price_data_vector.clear();
    
    for(int j = 0; j < RSI_data_vector.size(); j++){
        for(int k = 0; k < RSI_data_vector[j].size(); k++){
            delete[] RSI_data[j][k];
            RSI_data_vector[j][k].clear();
        }
        delete[] RSI_data[j];
        RSI_data_vector[j].clear();
    }

    delete[] RSI_data;
    RSI_data_vector.clear();
    delete[] data_copy;
}

string getOutputFilePath(string company_name, Date current_date, string mode, string file_dir, string type){
    return file_dir + "/" + company_name + "/" + type + "/" + mode + "/" + mode + "_" + current_date.getYear() + "_" + current_date.getMon() + ".csv";
}

void startTrain(RSIParticle &result, string company_name, CompanyData &companyData, int range_day_number){
    
    double *beta_ = new double[BIT_SIZE];
    
    RSIParticle expBest(range_day_number, BIT_SIZE, FUNDS, companyData);
    RSIParticle gBest(range_day_number, BIT_SIZE, FUNDS, companyData);
    RSIParticle gWorst(range_day_number, BIT_SIZE, FUNDS, companyData);
    RSIParticle pBest(range_day_number, BIT_SIZE, FUNDS, companyData);
    RSIParticle pWorst(range_day_number, BIT_SIZE, FUNDS, companyData);
    RSIParticle* rsi_particle_list = new RSIParticle[PARTICLENUMBER];
    initRSIParticle(rsi_particle_list, PARTICLENUMBER, range_day_number, companyData);
    for(int n = 0; n < EXPNUMBER; n++){
        
        cout << "___" << company_name << " : " << n << "___" << endl;
        gBest.init();
        gWorst.init();
        gBest.return_rate = -10000;
        gWorst.return_rate = DBL_MAX;
        initial(beta_, BIT_SIZE);
        for(int i = 0; i < ITERNUMBER; i++){
            pBest.init();
            pWorst.init();
            initRISParticle(rsi_particle_list);
            genParticle(rsi_particle_list, n, i, beta_);
            startTrade(rsi_particle_list, range_day_number, PARTICLENUMBER);
            recordGAnswer(rsi_particle_list, gBest, gWorst, pBest, pWorst);
            adjBeta(gBest, pWorst, beta_);
        }
        recordExpAnswer(expBest, gBest);
        
    }
    delete[] rsi_particle_list;
    expBest.print();
    delete[] beta_;
    result.copyP(expBest);
}

void startExhaustive(RSIParticle &result, string company_name, CompanyData &companyData, int range_day_number){
    
    RSIParticle expBest(range_day_number, BIT_SIZE, FUNDS, companyData);
    expBest.return_rate = -10000;
    RSIParticle* rsi_particle_list = new RSIParticle[1];
    initRSIParticle(rsi_particle_list, 1, range_day_number, companyData);
    int* temp_data = new int[BIT_SIZE];
    int EXHAUSTIVENUMBER = pow(2, BIT_SIZE);
    for(int n = 0; n < EXHAUSTIVENUMBER; n++){
        if(n % 100000 == 0){
            cout << "___" << company_name << " : " << n << "___" << endl;
        }
        rsi_particle_list[0].init();
       
        bool add = true;
        for(int j = 0; j < BIT_SIZE; j++){
            if (n == 0){
                temp_data[j] = 0;
            }
            if(temp_data[j] == 0 && add){
                add = false;
                temp_data[j] = 1;
                break;
            }else if(temp_data[j] == 1 && add){
                temp_data[j] = 0;
            }
        }
        
        for(int j = 0; j < BIT_SIZE; j++){
            rsi_particle_list[0].data[j] = temp_data[j];
        }
        rsi_particle_list[0].bitToDec();
        
        startTrade(rsi_particle_list, range_day_number, 1);
        recordExpAnswer(expBest, rsi_particle_list[0]);
    }
    delete[] rsi_particle_list;
    expBest.print();
    result.copyP(expBest);
}

int main(int argc, const char * argv[]) {
    
    
    vector<string> company_list = genCompanyList(COMPANY_PRICE_DIR);
    for(int c = 0; c < company_list.size(); c++){
        int day_number;
        string** price_data;
        string*** RSI_data;
        vector<vector<string>> price_data_vector;
        vector<vector<vector<string>>> RSI_data_vector;
        
        string temp = COMPANY_PRICE_DIR + "/" + company_list[c] + ".csv";
        readData(temp, price_data_vector);
        temp = RSI_DIR + "/" + company_list[c] + "/";
        readRSIData(temp, RSI_data_vector);
        
        price_data = vectorToArray(price_data_vector);
        RSI_data = vectorToArray(RSI_data_vector);
        day_number = price_data_vector.size() - 1;
        
        string *data_copy = new string[day_number];
        copyData(data_copy, price_data, day_number);
        
        for(int s = 4; s >= 0; s--){
            
            srand(343);
            double START, END;
            START = clock();
            
            Date current_date;
            Date finish_date;
            string TYPE;
            preSet(MODE, current_date, finish_date, s, TYPE);
            createDir(FILE_DIR, company_list[c], TYPE);
            
            temp = FILE_DIR + "/" + company_list[c] + "/" + TYPE + "/" + "total_data_" + MODE + ".csv";
            ofstream outfile_data;
            outfile_data.open(temp, ios::out);
            outfile_data << "Date,EXP,GEN,RSI number,Buy point,Sell point,Trade times,Return rate," << endl;
            
            do{
                int range_day_number;
                string start_date;
                string end_date;
                int start_index;
                int end_index;
                
                setWindow(MODE, start_date, end_date, start_index, end_index, current_date, finish_date, data_copy, price_data, day_number, range_day_number);
                CompanyData companyData;
                companyData.init(range_day_number);
                createCompanyData(companyData, price_data, RSI_data, start_index, end_index);
                range_day_number = companyData.day_number;
                
                cout << "______" << TYPE << " : " << start_date << " - " << end_date << "______" << endl;
                RSIParticle result(range_day_number, BIT_SIZE, FUNDS, companyData);
                
                if(MODE == "train"){
                    startTrain(result, company_list[c], companyData, range_day_number);
                }else if(MODE == "exhaustive"){
                    startExhaustive(result, company_list[c], companyData, range_day_number);
                }
                
                genTradeRecord(result, range_day_number);
                outputFile(result, getOutputFilePath(company_list[c], current_date, MODE, FILE_DIR, TYPE));
                recordData(result, outfile_data);
                current_date.slide();
            }while(finish_date >= current_date);
            END = clock();
            temp = FILE_DIR + "/" + company_list[c] + "/" + TYPE + "/" + "time_" + MODE + ".txt";
            recordCPUTime(START, END, temp);
            outfile_data.close();
        }
        releaseData(price_data_vector, RSI_data_vector, price_data, RSI_data, data_copy);
    }
    return 0;
}
