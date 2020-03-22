#ifndef FILE_IO_H
#define FILE_IO_H

#include <string>
#include <unordered_map>
#include "mesh.h"

class HashPair {
public:
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2> &p) const {
        auto hash1 = std::hash<T1>()(p.first),
                hash2 = std::hash<T2>()(p.second);
        return hash1 ^ hash2;
    }
}; 

class FileReader {
public:
    FileReader(const std::string &filePath);
    
    const std::unordered_map<int, VertexPtr> & vertices() const {
        return m_verts;
    }
    
    const std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> & edges() const {
        return m_edges;
    }
    
    const std::string & instanceName() const {
        return m_instanceName;
    }
    
private:
    std::unordered_map<int, VertexPtr> m_verts;
    std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> m_edges;
    std::string m_instanceName;
};

class FileWriter {
public:
    FileWriter(const std::string &filePath, const std::string &instanceName,
               const std::unordered_map<int, VertexPtr> &verts,
               const std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> &edges);
    
    void operator()() const;
    
private:
    std::string m_filePath;
    std::string m_instanceName;
    std::unordered_map<int, VertexPtr> m_verts;
    std::unordered_map<std::pair<int, int>, EdgePtr, HashPair> m_edges;
};

#endif /* FILE_IO_H */
