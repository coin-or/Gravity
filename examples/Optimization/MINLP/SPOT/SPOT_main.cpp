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

#ifdef USE_OUSTER
#include "helpers.h"
#include "ouster/impl/build.h"
#include "ouster/client.h"
#include "ouster/lidar_scan.h"
#include "ouster/os_pcap.h"
#include "ouster/types.h"
using namespace ouster;
img_t<double> get_x_in_image_form(const LidarScan& scan, bool destaggered,
                                  const sensor::sensor_info& info) {
    // For convenience, save w and h to variables
    const size_t w = info.format.columns_per_frame;
    const size_t h = info.format.pixels_per_column;

    // Get the XYZ in ouster::Points (n x 3 Eigen array) form
    XYZLut lut = make_xyz_lut(info);
    auto cloud = cartesian(scan.field(sensor::ChanField::RANGE), lut);

    // Access x and reshape as needed
    // Note that the values in cloud.col(0) are ordered
    auto x = Eigen::Map<const img_t<double>>(cloud.col(0).data(), h, w);
    auto x_destaggered = destagger<double>(x, info.format.pixel_shift_by_row);

    // Apply destagger if desired
    if (!destaggered) return x;
    return x_destaggered;
}

#endif

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
//        auto xyz = ipLayer->
    }
    
    // create the stats object
    pcpp::IPcapDevice::PcapStats stats;
    // read stats from reader and print them
    reader.getStatistics(stats);
    std::cout << "Read " << stats.packetsRecv << " packets successfully and " << stats.packetsDrop << " packets could not be read" << std::endl;
    // close the file
    reader.close();
#ifdef USE_OUSTER
    using namespace ouster::sensor;
    if (argc != 3) {
        std::cerr << "Version: " << ouster::SDK_VERSION_FULL << " ("
        << ouster::BUILD_SYSTEM << ")"
        << "\n\nUsage: lidar_scan_example <pcap_file> <json_file>"
        << std::endl;
        
        return argc == 1 ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    
    const std::string pcap_file = argv[1];
        const std::string json_file = argv[2];

        auto handle = sensor_utils::replay_initialize(pcap_file);
        auto info = sensor::metadata_from_json(json_file);

        size_t w = info.format.columns_per_frame;
        size_t h = info.format.pixels_per_column;

        auto scan = LidarScan(w, h, info.format.udp_profile_lidar);

        std::cerr << "Reading in scan from pcap..." << std::endl;
        get_complete_scan(handle, scan, info);

        // 1. Getting XYZ
        std::cerr << "1. Calculating 3d Points... " << std::endl;
        //! [doc-stag-cpp-xyz]
        XYZLut lut = make_xyz_lut(info);
        auto range = scan.field(sensor::ChanField::RANGE);
        auto cloud = cartesian(range, lut);
        //! [doc-etag-cpp-xyz]
        //
    std::cerr << "\nLet's see what's in this cloud...  " << endl;
    for(auto i = 0; i < cloud.size()-1; i++)
        DebugOn("(" << cloud(i, 0) << ", " << cloud(i, 1) << ", "
                  << cloud(i, 2) << ")" << std::endl);
        

        // 3. Destaggering
        // Fields come in w x h arrays, but they are staggered, so that a column
        // reflects the timestamp. To get each column to make visual sense,
        // destagger the image
        std::cerr
            << "\n3. Getting staggered and destaggered images of Reflectivity..."
            << std::endl;

        Eigen::Array<uint32_t, -1, -1, Eigen::RowMajor> reflectivity;

        if (info.format.udp_profile_lidar ==
            sensor::UDPProfileLidar::PROFILE_LIDAR_LEGACY) {
            reflectivity = scan.field(sensor::ChanField::REFLECTIVITY);
        } else if (info.format.udp_profile_lidar ==
                   sensor::UDPProfileLidar::PROFILE_RNG19_RFL8_SIG16_NIR16_DUAL) {
            reflectivity = scan.field<uint8_t>(sensor::ChanField::REFLECTIVITY)
                               .cast<uint32_t>();
        } else {  // legacy or single return profile
            reflectivity = scan.field<uint16_t>(sensor::ChanField::REFLECTIVITY)
                               .cast<uint32_t>();
        }

        //! [doc-stag-cpp-destagger]
        auto reflectivity_destaggered =
            destagger<uint32_t>(reflectivity, info.format.pixel_shift_by_row);
        //! [doc-etag-cpp-destagger]

        // 4. You can get XYZ in w x h arrays too
        std::cerr
            << "4. Getting staggered and destaggered images of X Coordinate..."
            << std::endl;
        auto x_image_staggered = get_x_in_image_form(scan, false, info);
        auto x_image_destaggered = get_x_in_image_form(scan, true, info);

        const auto print_row = std::min<size_t>(123, h - 3);
        const auto print_column = std::min<size_t>(1507, w / 2 + 5);

        const std::string point_string = "(" + std::to_string(print_row) + ", " +
                                         std::to_string(print_column) + ")";

        std::cerr << "In the staggered image, the point at " << point_string
                  << " has reflectivity " << reflectivity(print_row, print_row)
                  << " and an x coordinate of "
                  << x_image_staggered(print_column, print_column) << "."
                  << std::endl;
        std::cerr << "In the destagged image, the point at " << point_string
                  << " has reflectivity "
                  << reflectivity_destaggered(print_row, print_column)
                  << " and an x coordinate of "
                  << x_image_destaggered(print_row, print_column) << "."
                  << std::endl;
#endif
#else
    DebugOn("Gravity compiled without PcapPlusPlus!" << endl);
#endif
    
    return 0;
}
