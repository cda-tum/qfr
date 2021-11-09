/*
 * This file is part of JKQ QFR library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#ifndef QFR_Ecc_HPP
#define QFR_Ecc_HPP

#include "QuantumComputation.hpp"
#include "EccStatistics.hpp"

class Ecc {
public:

    enum ID {
		Id, Q3Shor, Q9Shor, Q7Steane, Q9Surface, QxCustom
	};
    struct Info {
        ID id;
        int nRedundantQubits;
        int nClassicalBitsPerQubit;
        std::string name;
    };

    const Info ecc;

	Ecc(Info ecc, qc::QuantumComputation& qc, int measureFrequency);
	virtual ~Ecc() = default;

	qc::QuantumComputation& apply();

    virtual std::ostream& printResult(std::ostream& out);

    virtual void dumpResult(const std::string& outputFilename);

	virtual void dumpResult(const std::string& outputFilename, qc::Format format) {
		size_t slash = outputFilename.find_last_of('/');
		size_t dot = outputFilename.find_last_of('.');
		statistics.outputName = outputFilename.substr(slash+1, dot-slash-1);
		qcMapped.dump(outputFilename, format);
	}

	virtual void dumpResult(std::ostream& os, qc::Format format) {
		qcMapped.dump(os, format);
	}


protected:

    qc::QuantumComputation& qc;
	qc::QuantumComputation qcMapped;
	EccStatistics statistics{};
	const int measureFrequency;
	bool decodingDone;

	virtual void writeEncoding()=0;

	virtual void measureAndCorrect()=0;

	virtual void writeDecoding()=0;

	virtual void mapGate(const std::unique_ptr<qc::Operation> &gate)=0;

	void gateNotAvailableError(const std::unique_ptr<qc::Operation> &gate);

	void writeToffoli(int target, int c1, bool p1, int c2, bool p2);
	void writeZToffoli(int target, int c1, bool p1, int c2, bool p2);

};


#endif //QFR_Ecc_HPP
