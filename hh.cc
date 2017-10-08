
/********************************************************************
Approximate frequent items in a data stream
G. Cormode 2002, 2003,2005
R. Ben Basat and G. Einziger, 2016-17

Last modified: 2017-10
*********************************************************************/

//#define PCAP
// you need to include libraries wpcap.lib ws2_32.lib to compile using PCAP.

#include "prng.h"
#include "lossycount.h"
#include "countmin.h"
#include "FAST.hpp"
#include <fstream>
#include <math.h>       /* sqrt */
#ifdef __GNUC__
#include <sys/time.h>
#include <time.h>
#endif
#define RUN_FAST 1 // Run FAST with gamma=4 (or as supplied by a command line argument)
#define RUN_FAST2 0 // FAST with gamma = 0.25
#define RUN_CMS 0 // Run the Count Min Sketch
#define RUN_CMH 0 // Run the Heirarchical Count Min Sketch
#define RUN_CS 0 // Run the Count Sketch
#define RUN_SSH 0 // Run the heap-based implementation of Space Saving
#define RUN_SSL 0 // Run the linked-lists-based implementation of Space Saving

#define RUN_MSE 0 // Compute the Root Mean Square Error (RMSE) of the algorithm instead of runtime
#define MSE_OF_CS 0
#define PROFILING 0


#define SAME_NR_COUNTERS 0
/******************************************************************/

#ifdef PCAP
#include <pcap.h>

#ifdef __GNUC__
#include <arpa/inet.h>
#endif

typedef struct ip_address{
    u_char byte1;
    u_char byte2;
    u_char byte3;
    u_char byte4;
}ip_address;

/* IPv4 header */
typedef struct ip_header{
    u_char  ver_ihl;        // Version (4 bits) + Internet header length (4 bits)
    u_char  tos;            // Type of service 
    u_short tlen;           // Total length 
    u_short identification; // Identification
    u_short flags_fo;       // Flags (3 bits) + Fragment offset (13 bits)
    u_char  ttl;            // Time to live
    u_char  proto;          // Protocol
    u_short crc;            // Header checksum
    ip_address  saddr;      // Source address
    ip_address  daddr;      // Destination address
    u_int   op_pad;         // Option + Padding
}ip_header;

/* UDP header*/
typedef struct udp_header{
    u_short sport;          // Source port
    u_short dport;          // Destination port
    u_short len;            // Datagram length
    u_short crc;            // Checksum
}udp_header;

typedef struct tcp_header{
	u_short sport;
	u_short dport;
}tcp_header;
#else
#ifdef _MSC_VER
#include <windows.h>
#endif
#endif

class Stats
{
public:
	Stats() : dU(0.0), dQ(0.0), dP(0.0), dR(0.0), dF(0.0), dF2(0.0) {}

	double dU, dQ;
	double dUMax, dUMin;
	double dP, dR, dF, dF2;
	std::multiset<double> P, R, F, F2;
};

void usage()
{
	std::cerr
		<< "Usage:\n"
		<< "  -np		number of packets\n"
		<< "  -r		number of runs\n"
		<< "  -phi		phi\n"
		<< "  -d		depth\n"
		<< "  -g		granularity\n"
		<< "  -gamma    DIM-SUM coefficient\n"
#ifdef PCAP
		<< "  -pf   pcap file name\n"
		<< "  -f    pcap filter\n"
#else
		<< "  -z    skew\n"
#endif
		<< std::endl;
}

void StartTheClock(uint64_t& s)
{
#ifdef _MSC_VER
	FILETIME ft;
    LARGE_INTEGER li;

	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;
	s = (uint64_t) (li.QuadPart / 10000);
#else
	struct timeval tv;
	gettimeofday(&tv, 0);
	s = (1000 * tv.tv_sec) + (tv.tv_usec / 1000);
#endif
}

// returns milliseconds.
uint64_t StopTheClock(uint64_t s)
{
#ifdef _MSC_VER
    FILETIME ft;
    LARGE_INTEGER li;
    uint64_t t;

	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;
	t = (uint64_t) (li.QuadPart / 10000);
	return t - s;
#else
	struct timeval tv;
	gettimeofday(&tv, 0);
	return (1000 * tv.tv_sec) + (tv.tv_usec / 1000) - s;
#endif
}

void CheckOutput(std::map<uint32_t, uint32_t>& res, uint64_t thresh, size_t hh, Stats& S, const std::vector<uint32_t>& exact)
{
	return;
	if (res.empty())
	{
		S.F.insert(0.0);
		S.F2.insert(0.0);
		S.P.insert(100.0);
		S.dP += 100.0;

		if (hh == 0)
		{
			S.R.insert(100.0);
			S.dR += 100.0;
		}
		else
			S.R.insert(0.0);

		return;
	}

	size_t correct = 0;
	size_t claimed = res.size();
	size_t falsepositives = 0;
	double e = 0.0, e2 = 0.0;

	std::map<uint32_t, uint32_t>::iterator it;
	for (it = res.begin(); it != res.end(); ++it)
	{
		if (exact[it->first] >= thresh)
		{
			++correct;
			uint32_t ex = exact[it->first];
			double diff = (ex > it->second) ? ex - it->second : it->second - ex;
			e += diff / ex;
		}
		else
		{
			++falsepositives;
			uint32_t ex = exact[it->first];
			double diff = (ex > it->second) ? ex - it->second : it->second - ex;
			e2 += diff / ex;
		}
	}

	if (correct != 0)
	{
		e /= correct;
		S.F.insert(e);
		S.dF += e;
	}
	else
		S.F.insert(0.0);

	if (falsepositives != 0)
	{
		e2 /= falsepositives;
		S.F2.insert(e2);
		S.dF2 += e2;
	}
	else
		S.F2.insert(0.0);

	double r = 100.0;
	if (hh != 0) r = 100.0 *((double) correct) / ((double) hh);

	double p = 100.0 *((double) correct) / ((double) claimed);

	S.R.insert(r);
	S.dR += r;
	S.P.insert(p);
	S.dP += p;
}

void PrintTimes(char* title, std::vector<uint64_t> times) {
	std::cout << title;
	for (auto const& t : times) {
		std::cout << "\t" << t;
	}
	std::cout << std::endl;
}

void PrintOutput(char* title, size_t size, const Stats& S, size_t u32NumberOfPackets, int reps, double RMSE)
{
	double p5th = -1.0, p95th = -1.0, r5th = -1.0, r95th = -1.0, f5th = -1.0, f95th = -1.0, f25th = -1.0, f295th = -1.0;
	size_t i5, i95;
	std::multiset<double>::const_iterator it;

	if (! S.P.empty())
	{
		it = S.P.begin();
		i5 = (size_t) (S.P.size() * 0.05);
		for (size_t i = 0; i < i5; ++i) ++it;
		p5th = *it;
		i95 = (size_t) (S.P.size() * 0.95);
		for (size_t i = 0; i < (i95 - i5); ++i) ++it;
		p95th = *it;
	}

	if (! S.R.empty())
	{
		it = S.R.begin();
		i5 = S.R.size() * 0.05;
		for (size_t i = 0; i < i5; ++i) ++it;
		r5th = *it;
		i95 = S.R.size() * 0.95;
		for (size_t i = 0; i < (i95 - i5); ++i) ++it;
		r95th = *it;
	}

	if (! S.F.empty())
	{
		it = S.F.begin();
		i5 = S.F.size() * 0.05;
		for (size_t i = 0; i < i5; ++i) ++it;
		f5th = *it;
		i95 = S.F.size() * 0.95;
		for (size_t i = 0; i < (i95 - i5); ++i) ++it;
		f95th = *it;
	}

	if (! S.F2.empty())
	{
		it = S.F2.begin();
		i5 = S.F2.size() * 0.05;
		for (size_t i = 0; i < i5; ++i) ++it;
		f25th = *it;
		i95 = S.F2.size() * 0.95;
		for (size_t i = 0; i < (i95 - i5); ++i) ++it;
		f295th = *it;
	}
#if RUN_MSE
	printf("%s\t%1.2f\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n",
		title, RMSE , size,
		S.dR / S.R.size(), r5th, r95th,
		S.dP / S.P.size(), p5th, p95th,
		S.dF / S.F.size(), f5th, f95th,
		S.dF2 / S.F2.size(), f25th, f295th);
#else
	printf("%s\t%1.2f\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n",
		title, u32NumberOfPackets  * reps / S.dU, size,
		S.dR / S.R.size(), r5th, r95th,
		S.dP / S.P.size(), p5th, p95th,
		S.dF / S.F.size(), f5th, f95th,
		S.dF2 / S.F2.size(), f25th, f295th
	);
#endif

}

size_t RunExact(uint64_t thresh, std::vector<uint32_t>& exact)
{
	size_t hh = 0;

	for (size_t i = 0; i < exact.size(); ++i)
		if (exact[i] >= thresh) ++hh;

	return hh;
}

/******************************************************************/

int main(int argc, char **argv) 
{
	size_t stNumberOfPackets = 100000000;
	size_t stRuns = 20;
	double dPhi = 1. / (1 << 20); // The error parameter. Allowed error = M*N*dPhi
	double gamma = 4.; // The gamma using which FAST is initialized
	int M = 100; // A bound on the maximal packet size
	int reps = 1; // How many times to repeat the measurement
	uint32_t u32Depth = 10; // How many rows in Count Min / Count Sketch
	uint32_t u32Granularity = 8;
	std::string file = ""; // The input trace file
	bool timeLaspe = false;
	bool exactIsValid = true;
	double secondGamma = 0.25;
	bool RunMSE = false;
	double MSE = 0;
	double MSE_Error = 0;
#ifdef PCAP
	std::string sFilter = "ip and (tcp or udp)";
	std::string sPcapFile = "header_gs10";
#else
	double dSkew = 1.0;
#endif
	//Parsing command arguments
	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "-np") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing number of packets." << std::endl;
				return -1;
			}
			stNumberOfPackets = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-r") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing number of runs." << std::endl;
				return -1;
			}
			stRuns = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-d") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing depth." << std::endl;
				return -1;
			}
			u32Depth = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-g") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing granularity." << std::endl;
				return -1;
			}
			u32Granularity = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-phi") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing phi." << std::endl;
				return -1;
			}
			dPhi = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-skew") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing skew." << std::endl;
				return -1;
			}
			dSkew = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-f") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing file name." << std::endl;
				return -1;
			}
			file = std::string(argv[i]);
			std::cerr << "File name =" << file << std::endl;
		}
		else if (strcmp(argv[i], "-t") == 0)
		{
			timeLaspe = true;
		}
		else if (strcmp(argv[i], "-gamma") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing gamma." << std::endl;
				return -1;
			}
			gamma = atof(argv[i]);
			std::cerr << "gamma =" << gamma << std::endl;
		}
		else if (strcmp(argv[i], "-M") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing M." << std::endl;
				return -1;
			}
			M = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-reps") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing reps." << std::endl;
				return -1;
			}
			reps = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-MSE") == 0)
		{
			i++;
			RunMSE = true;
		}
#ifdef PCAP
		else if (strcmp(argv[i], "-pf") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing pcap file name." << std::endl;
				return -1;
			}
			sPcapFile = std::string(argv[i]);
		}
		else if (strcmp(argv[i], "-f") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing pcap filter." << std::endl;
				return -1;
			}
			sFilter = std::string(argv[i]);
		}
#else
		else if (strcmp(argv[i], "-z") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cerr << "Missing skew parameter." << std::endl;
				return -1;
			}
			dSkew = atof(argv[i]);
		}
#endif
		else
		{
			usage();
			return -1;
		}
	}

#if SAME_NR_COUNTERS==1 
	uint32_t u32Width = ceil( (1. / dPhi ) / u32Depth);
#else
	uint32_t u32Width = ceil(2.7183 / dPhi);
#endif

	prng_type * prng;
	prng=prng_Init(44545,2);
	int64_t a = (int64_t) (prng_int(prng)% MOD);
	int64_t b = (int64_t) (prng_int(prng)% MOD);
	prng_Destroy(prng);

	uint32_t u32DomainSize = 1048575;
	std::vector<uint32_t> exact(u32DomainSize + 1, 0), exactMSE(u32DomainSize + 1, 0);
#ifdef PCAP
	//Open the capture file
	pcap_t *fp;
	char errbuf[PCAP_ERRBUF_SIZE];
	if ((fp = pcap_open_offline(sPcapFile.c_str(), errbuf)) == 0)
	{
		std::cerr << "Unable to open file." << std::endl;
		exit(1);
	}

    //compile the filter
	struct bpf_program fcode;
	if (pcap_compile(fp, &fcode, const_cast<char*>(sFilter.c_str()), 1, 0xffffff) < 0 )
    {
		std::cerr << "Unable to compile the packet filter. Check the syntax" << std::endl;
        exit(1);
    }

    //set the filter
    if (pcap_setfilter(fp, &fcode) < 0)
    {
		std::cerr << "Error setting the filter" << std::endl;
        exit(1);
    }

	struct pcap_pkthdr *header;
	const u_char *pkt_data;
	int res;
#endif
	Stats SLS, SCMH, SCM, SCMS, SCCFC, SALS, SLCL, SLCU, SFAST, SFAST2;
	std::vector<uint64_t> TLS, TCMH, TCM, TCCFC, TALS, TLCL,TLCU;
	CM_type* cm = NULL;
	CMH_type* cmh = NULL;
	CS_type *cs = NULL;
	LCL_type* lcl = NULL;
	LCU_type* lcu = NULL;

	int MSE_Packets = 0;
	bool MSEExactIsValid = true;
	if (RUN_MSE && RUN_FAST)
		throw std::exception();

	FAST *fast, *fast2, *FASTMSE;
	if (RUN_FAST)
		fast = new FAST(dPhi, M, gamma);
	if (RUN_FAST2)
		fast2 = new FAST(dPhi, M, 0.25f);
	clock_t begin, end;
	if (RUN_SSH)
		lcl = LCL_Init(dPhi);
	if (RUN_SSL)
		lcu = LCU_Init(dPhi);
	if (RUN_CMS)
		cm = CM_Init(u32Width, u32Depth, 32);
	if (RUN_CMH)
		cmh = CMH_Init(u32Width, u32Depth, 32, u32Granularity);
	if (RUN_CS)
		cs = CS_Init(u32Width*u32Width, u32Depth, 32);

	if (RUN_MSE)
		FASTMSE = new FAST(dPhi, M, gamma);
		
	std::vector<uint32_t> data;
	std::vector<uint32_t> values;
#ifndef PCAP
	Tools::Random r = Tools::Random(0xF4A54B);
	Tools::PRGZipf zipf = Tools::PRGZipf(0, u32DomainSize, dSkew, &r);
#endif

	size_t stCount = 0;
#ifdef PCAP
	while((res = pcap_next_ex( fp, &header, &pkt_data)) >= 0 && stCount < stNumberOfPackets)
#else
	if (file != "") {
		uint64_t total = 0;
		std::ifstream f(file);
		int id, length;
		while (f >> id >> length) {
			assert(length > 0);
			if (total < 0) {
				std::cerr << "Why is total negative? " << total<<std::endl;
				break;
			}
			
			if (total >= 2000000) {
				stNumberOfPackets = total; 
				break; 
			}
				
			data.push_back(id);
			values.push_back(length);
			total ++;
		}
		if (data.size() / stRuns < 400000)
			stRuns = data.size() / 400000;
		std::cerr << "Finished loading file. Total number of bytes: " << total << std::endl;
	}
	else {
		for (int i = 0; i < stNumberOfPackets; ++i)
#endif
		{
			++stCount;
			if (stCount % 5000000 == 0)
				std::cerr << stCount << std::endl;
#ifdef PCAP
			ip_header *ih;
			udp_header *uh;
			tcp_header *th;
			u_int ip_len;
			u_short sport, dport;

			//retireve the position of the ip header
			ih = (ip_header *)(pkt_data + 14); //length of ethernet header
			ip_len = (ih->ver_ihl & 0xf) * 4;

			if (ih->proto == 6)
			{
				th = (tcp_header *)((u_char*)ih + ip_len);
				sport = ntohs(th->sport);
				dport = ntohs(th->dport);
			}
			else if (ih->proto == 17)
			{
				uh = (udp_header *)((u_char*)ih + ip_len);
				sport = ntohs(uh->sport);
				dport = ntohs(uh->dport);
			}

			uint32_t v;
			memcpy(&v, &(ih->daddr), 4);
#else
			uint32_t v = zipf.nextLong();
#endif
			uint32_t value = v & u32DomainSize;
			if (value > 0) {
				data.push_back(value);
				values.push_back(1);
			}
			else {
				data.push_back(-value);
				values.push_back(1);
			}
		}
	}
	size_t stRunSize = data.size() / stRuns;
	stNumberOfPackets = data.size();
	for (size_t iteration = 0; iteration < reps; ++iteration)
	{
		size_t stStreamPos = 0;
		uint64_t nsecs;
		uint64_t t;
		unsigned long long total = 0;
		for (size_t run = 1; run <= stRuns; ++run) // stRuns
		{
			bool stop = false;
			if (exactIsValid) {
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					assert(values[i] > 0);
					total += abs((int)values[i]);
					if (total >= 0x7FFFFFFF) {
						exactIsValid = false;
						if (!PROFILING)
							std::cerr << "Error! Total number of bytes is " << total << std::endl;
						break;
					}
					exact[data[i]] += values[i];
					if (exact[data[i]] > 0x7FFFFFFF) {
						exactIsValid = false;
						stop = true;
						std::cerr << "Strange. Value is too large " << exact[data[i]] << " after addding " << values[i] << std::endl;
						break;
					}
				}
			}
			if (stop) {
				break;
			}

			if (RUN_CMS) {
				begin = clock();
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					CM_Update(cm, data[i], values[i]);
				}
				end = clock();
				SCM.dU += t = end - begin;
				TCM.push_back(t);
				
			}
			if (RUN_CMH) {
				begin = clock();
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					CMH_Update(cmh, data[i], values[i]);
				}
				end = clock();
				SCMH.dU += t = end - begin;
				TCMH.push_back(t);
			}
			if (RUN_CS) {
				begin = clock();
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					CS_Update(cs, data[i], values[i]);
				}
				end = clock();
				SCCFC.dU += t = end - begin;
				TCCFC.push_back(t);
			}
			if (RUN_SSH) {
				begin = clock();
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					LCL_Update(lcl, data[i], values[i]);
				}
				end = clock();
				SLCL.dU += t = end - begin;
				TLCL.push_back(t);
			}
			if (RUN_SSL) {
				begin = clock();
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					LCU_Update(lcu, data[i]);
				}
				end = clock();
				SLCU.dU += t = end - begin;
				TLCU.push_back(t);
			}
			
			if (RUN_FAST){
				begin = clock();
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					fast->update(data[i], values[i]);
				}
				end = clock();
				SFAST.dU += end - begin;
			}
#if RUN_MSE
			for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
			{
				if (MSEExactIsValid){
#if MSE_OF_CS==1
					long err = CS_Query(cs, data[i]) - exactMSE[data[i]];
#else
					long err = FASTMSE->query(data[i]) - exactMSE[data[i]];
#endif
					double prevMSE = MSE;
					MSE_Error += (double) err * err;
					MSE += MSE_Error;
					MSE_Error = MSE_Error - (MSE - prevMSE);
					exactMSE[data[i]] += values[i];
					if (exactMSE[data[i]] > 0x7FFFFFFF)
						MSEExactIsValid = false;
					MSE_Packets += 1;
					if (MSEExactIsValid) {
						FASTMSE->update(data[i], values[i]);
						assert(FASTMSE->query(data[i]) > 0);
#if  _DEBUG
						if (UNDER_ESTIMATOR) {
							if (exactMSE[data[i]] < FASTMSE->query(data[i]))
								std::cerr << data[i] << ", " << values[i] << ", " << exactMSE[data[i]] << ", Error = " << err << "," << i << std::endl;
							assert(exactMSE[data[i]] >= FASTMSE->query(data[i]));
							if (exactMSE[data[i]] > FASTMSE->query(data[i]) + i*dPhi*M)
								std::cerr << data[i] << ", " << values[i] << ", " << exactMSE[data[i]] << ", Error = " << err << "," << i << std::endl;
							assert(exactMSE[data[i]] <= FASTMSE->query(data[i]) + i*dPhi*M);
						}
						else {
							assert(exactMSE[data[i]] <= FASTMSE->query(data[i]));
							assert(exactMSE[data[i]] >= FASTMSE->query(data[i]) - i*dPhi*M);
						}
#endif
					}
				}
			}
#endif
			if (RUN_FAST2) {
				begin = clock();
				for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
				{
					fast2->update(data[i], values[i]);
				}
				end = clock();
				SFAST2.dU += end - begin;
			}
			uint64_t thresh = static_cast<uint64_t>(floor(dPhi*total) + 1);
			size_t hh = RunExact(thresh, exact);
			if (!PROFILING)
				std::cerr << "Run: " << run << ", Exact: " << hh << std::endl;

			std::map<uint32_t, uint32_t> res;
			stStreamPos += stRunSize;
		}
	}
	if (timeLaspe) {
		PrintTimes("CM", TCMH);
		PrintTimes("CS", TCCFC);
		PrintTimes("SSH", TLCL);
	}
	else {
		if (RUN_MSE)
			printf("\nMethod\tRMSE\tSpace\tRecall\t5th\t95th\tPrecis\t5th\t95th\tFreq RE\t5th\t95th\n");
		else
			printf("\nMethod\tUpdates/ms\tSpace\tRecall\t5th\t95th\tPrecis\t5th\t95th\tFreq RE\t5th\t95th\n");
		if (RUN_CMS)
			PrintOutput("CM", CM_Size(cm), SCM, stNumberOfPackets, reps, 0.0);
		if (RUN_CMH)
			PrintOutput("CMH", CMH_Size(cmh), SCMH, stNumberOfPackets, reps, 0.0);
		if (RUN_CS)
			PrintOutput("CCFC", CS_Size(cs), SCCFC, stNumberOfPackets, reps, 0.0);
		if (RUN_SSH)
			PrintOutput("SSH", LCL_Size(lcl), SLCL, stNumberOfPackets, reps, 0.0);
		if (RUN_SSL)
			PrintOutput("SSL", LCU_Size(lcu), SLCU, stNumberOfPackets, reps, 0.0);
		if (RUN_FAST)
			PrintOutput("FAST", fast->size(), SFAST, stNumberOfPackets, reps, 0.0);
		if (RUN_FAST2)
			PrintOutput("FAST2", fast2->size(), SFAST2, stNumberOfPackets, reps, 0.0);
		if (RUN_MSE)
			PrintOutput("FAST", fast->size(), SFAST, stNumberOfPackets, reps, sqrt(MSE / MSE_Packets));
	}
	
	if (RUN_CMS)
		CM_Destroy(cm);
	if (RUN_CS)
		CS_Destroy(cs);
	if (RUN_SSH)
		LCL_Destroy(lcl);	
	printf("\n");
	return 0;
}

