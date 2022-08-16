////
////
////  Gravity
////
////  Created by Hassan Hijazi on August 12 2022.
///
/// Optimization for the SPOT robot Lidar scanner
////
////
//
//

#include <iostream>
#include <gravity/solver.h>
#ifdef USE_PcapPlusPlus
#include "IPv4Layer.h"
#include "Packet.h"
#include "PcapFileDevice.h"
#include "SystemUtils.h"
#include "EthLayer.h"
#include "IPv4Layer.h"
#include "TcpLayer.h"
#include "HttpLayer.h"



std::string getProtocolTypeAsString(pcpp::ProtocolType protocolType)
{
    switch (protocolType)
    {
        case pcpp::Ethernet:
            return "Ethernet";
        case pcpp::IPv4:
            return "IPv4";
        case pcpp::TCP:
            return "TCP";
        case pcpp::HTTPRequest:
        case pcpp::HTTPResponse:
            return "HTTP";
        default:
            return "Unknown";
    }
}

std::string printTcpFlags(pcpp::TcpLayer* tcpLayer)
{
    std::string result = "";
    if (tcpLayer->getTcpHeader()->synFlag == 1)
        result += "SYN ";
    if (tcpLayer->getTcpHeader()->ackFlag == 1)
        result += "ACK ";
    if (tcpLayer->getTcpHeader()->pshFlag == 1)
        result += "PSH ";
    if (tcpLayer->getTcpHeader()->cwrFlag == 1)
        result += "CWR ";
    if (tcpLayer->getTcpHeader()->urgFlag == 1)
        result += "URG ";
    if (tcpLayer->getTcpHeader()->eceFlag == 1)
        result += "ECE ";
    if (tcpLayer->getTcpHeader()->rstFlag == 1)
        result += "RST ";
    if (tcpLayer->getTcpHeader()->finFlag == 1)
        result += "FIN ";
    
    return result;
}

std::string printTcpOptionType(pcpp::TcpOptionType optionType)
{
    switch (optionType)
    {
        case pcpp::PCPP_TCPOPT_NOP:
            return "NOP";
        case pcpp::PCPP_TCPOPT_TIMESTAMP:
            return "Timestamp";
        default:
            return "Other";
    }
}

std::string printHttpMethod(pcpp::HttpRequestLayer::HttpMethod httpMethod)
{
    switch (httpMethod)
    {
        case pcpp::HttpRequestLayer::HttpGET:
            return "GET";
        case pcpp::HttpRequestLayer::HttpPOST:
            return "POST";
        default:
            return "Other";
    }
}
#endif

int main(int argc, char* argv[])
{
    int ethPacketCount;
    int ipv4PacketCount;
    int ipv6PacketCount;
    int tcpPacketCount;
    int udpPacketCount;
    int dnsPacketCount;
    int httpPacketCount;
    int sslPacketCount;
    
    string fname = string(prj_dir)+"/data_sets/SPOT/1_packet.pcap";
    if(argc==2){
        fname=argv[1];
    }
#ifdef USE_PcapPlusPlus
    // open a pcap file for reading
    pcpp::PcapFileReaderDevice reader(fname);
    if (!reader.open())
    {
        std::cerr << "Error opening the pcap file" << std::endl;
        return 1;
    }
    
    // read the first (and only) packet from the file
    pcpp::RawPacket rawPacket;
//    if (!reader.getNextPacket(rawPacket))
//    {
//        std::cerr << "Couldn't read the first packet in the file" << std::endl;
//        return 1;
//    }
    
    
    while (reader.getNextPacket(rawPacket))
    {
        // parse the raw packet into a parsed packet
        pcpp::Packet parsedPacket(&rawPacket);
        // let's get the IPv4 layer
        pcpp::IPv4Layer* ipLayer = parsedPacket.getLayerOfType<pcpp::IPv4Layer>();
        if (ipLayer == NULL)
        {
            std::cerr << "Something went wrong, couldn't find IPv4 layer" << std::endl;
            return 1;
        }
        
        // print source and dest IP addresses, IP ID and TTL
        std::cout << std::endl
        << "Source IP address: " << ipLayer->getSrcIPAddress() << std::endl
        << "Destination IP address: " << ipLayer->getDstIPAddress() << std::endl
        << "IP ID: 0x" << std::hex << pcpp::netToHost16(ipLayer->getIPv4Header()->ipId) << std::endl
        << "TTL: " << std::dec << (int)ipLayer->getIPv4Header()->timeToLive << std::endl;
        auto xyz = ipLayer->
    }
    
    // create the stats object
    pcpp::IPcapDevice::PcapStats stats;
    // read stats from reader and print them
    reader.getStatistics(stats);
    std::cout << "Read " << stats.packetsRecv << " packets successfully and " << stats.packetsDrop << " packets could not be read" << std::endl;
    // close the file
    reader.close();
#else
    DebugOn("Gravity compiled without PcapPlusPlus!" << endl);
#endif
    
    return 0;
}
