#include <fstream>
#include <iostream>
#include <regex>
#include <jsoncpp/json/value.h>
#include <jsoncpp/json/json.h>
#include <jsoncpp/json/writer.h>
#include "file_io.h"

FileReader::FileReader(const std::string &filePath) {
    std::cout << "Opening file " << filePath << " ...\n";
    std::ifstream ifstr(filePath, std::ifstream::binary);

    if (!ifstr) {
        std::cout << "Error:  Could not open input file.\n";
    } else {
        const auto suffix = filePath.substr(filePath.find_last_of(".") + 1);
        const auto fSplit
                = [] (const std::string &str, const std::regex & re) {
                    return std::vector<std::string>{std::sregex_token_iterator(str.begin(), str.end(), re, -1),
                                                    std::sregex_token_iterator()};
                };

        if (suffix == "json") {
            Json::Value root;
            ifstr >> root;
            const Json::Value &pnts = root["points"];
            m_instanceName = root["name"].asString();

            for (const auto &pnt : pnts) {
                int id = pnt["i"].asInt();
                long long x = pnt["x"].asInt64(), y = pnt["y"].asInt64();
                m_verts[id] = std::make_shared<Vertex>(id, GmPnt(x, y));
            }
        } else if (suffix == "pts") {
            int id = 0;
            std::string line;
            while (std::getline(ifstr, line)) {
                const auto strs = fSplit(line, std::regex("\\s+"));
                
                if (strs.size() == 2) {
                    double x = std::stod(strs.at(0)), y = std::stod(strs.at(1));
                    m_verts[id] = std::make_shared<Vertex>(id, GmPnt(x, y));
                    ++id;
                }
            }
        } 
        
        /*else if (suffix == "obj") {
            const auto fSplit =
                    [] (const std::string &str, const std::regex & re) {
                        return std::vector<std::string>{std::sregex_token_iterator(str.begin(), str.end(), re, -1),
                                                        std::sregex_token_iterator()};
                    };

            int id = 1;
            std::string line;
            while (std::getline(ifstr, line)) {
                const auto strs = fSplit(line, std::regex("\\s+"));
                const auto c = strs.at(0);

                if (c == "v") {
                    long long x = std::stoll(strs.at(1)), y = std::stoll(strs.at(2));
                    m_verts[id++] = std::make_shared<Vertex>(id, GmPnt(x, y));
                } else if (c == "f") {
                    if (strs.size() > 2) {
                        for (size_t i = 1; i < strs.size() - 1; i++) {
                            int from = std::stoi(strs.at(i)), 
                                    to = std::stoi(strs.at(i + 1));
                            const auto edgeId = std::minmax({from, to});

                            if (m_edges.find(edgeId) == m_edges.end()) {
                                const auto v = m_verts.at(from), w = m_verts.at(to);
                                m_edges[edgeId] = std::make_shared<Edge>(v, w, true);
                            } else {
                                m_edges[edgeId]->setIsOnConvexHull(false);
                            }
                        }
                    }
                }
            }
        }*/
    }
}

FileWriter::FileWriter(const std::string &filePath, const std::string &instanceName,
                       const std::unordered_map<int, VertexPtr> &verts,
                       const std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> &edges)
: m_filePath(filePath)
, m_instanceName(instanceName)
, m_verts(verts)
, m_edges(edges) {
}

void FileWriter::operator()() const {
    std::cout << "Writing results to " << m_filePath << " ...\n";
    std::ofstream ofstr(m_filePath);

    if (!ofstr) {
        std::cout << "Error:  Could not open input file.\n";
    } else {
        const auto suffix = m_filePath.substr(m_filePath.find_last_of(".") + 1);

        if (suffix == "json") {
            Json::Value decomp;
            decomp["type"] = "Solution";
            decomp["instance_name"] = m_instanceName;

            Json::Value edges(Json::arrayValue);
            for (const auto &entry : m_edges) {
                const auto edgeId = entry.first;
                Json::Value edge;
                edge["i"] = edgeId.first;
                edge["j"] = edgeId.second;
                edges.append(edge);
            }

            decomp["edges"] = edges;
            ofstr << decomp;
        } else if (suffix == "obj") {
            for (const auto &entry : m_verts) {
                const auto pnt = entry.second->pnt();
                ofstr << "v " << pnt.x() << " " << pnt.y() << "\n";
            }
            
            for (const auto &entry : m_edges) {
                const auto edgeId = entry.first;
                ofstr << "l " << (edgeId.first + 1) << " " 
                        << (edgeId.second + 1) << "\n";
            }
            
            ofstr << "\n";
        }
    }
}
