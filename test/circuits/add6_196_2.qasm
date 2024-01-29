OPENQASM 2.0;
include "qelib1.inc";
qreg q[19];
creg c[19];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
mcx q[16],q[10],q[9],q[2];
ccx q[16],q[10],q[4];
ccx q[16],q[10],q[3];
ccx q[17],q[11],q[5];
ccx q[17],q[11],q[4];
ccx q[17],q[11],q[3];
ccx q[15],q[9],q[3];
ccx q[15],q[9],q[2];
ccx q[15],q[9],q[1];
ccx q[15],q[9],q[0];
mcx q[15],q[9],q[8],q[1];
mcx q[15],q[9],q[8],q[0];
mcx q[14],q[8],q[7],q[0];
mcx q[16],q[15],q[10],q[2];
mcx q[16],q[15],q[10],q[1];
mcx q[16],q[15],q[14],q[10],q[1];
mcx q[18],q[12],q[11],q[4];
mcx q[18],q[12],q[11],q[3];
mcx q[18],q[12],q[11],q[2];
mcx q[18],q[12],q[11],q[1];
mcx q[17],q[11],q[10],q[3];
mcx q[14],q[13],q[8],q[0];
ccx q[14],q[8],q[2];
ccx q[14],q[8],q[1];
ccx q[13],q[7],q[1];
ccx q[13],q[7],q[0];
x q[16];
mcx q[17],q[16],q[15],q[11],q[8],q[1];
mcx q[17],q[16],q[11],q[9],q[8],q[1];
mcx q[17],q[16],q[11],q[3];
mcx q[18],q[16],q[12],q[11],q[9],q[2];
mcx q[18],q[16],q[12],q[11],q[9],q[1];
mcx q[18],q[17],q[16],q[12],q[9],q[2];
mcx q[18],q[17],q[16],q[12],q[9],q[1];
x q[10];
mcx q[18],q[15],q[14],q[12],q[11],q[10],q[7],q[0];
mcx q[18],q[17],q[15],q[14],q[12],q[10],q[7],q[0];
mcx q[17],q[11],q[10],q[9],q[8],q[1];
mcx q[17],q[15],q[11],q[10],q[8],q[1];
ccx q[16],q[10],q[4];
x q[17];
mcx q[18],q[17],q[12],q[4];
mcx q[18],q[17],q[15],q[14],q[12],q[10],q[1];
mcx q[18],q[17],q[15],q[12],q[10],q[8],q[7],q[0];
mcx q[18],q[17],q[12],q[10],q[9],q[2];
mcx q[18],q[17],q[12],q[10],q[9],q[1];
mcx q[18],q[17],q[12],q[10],q[9],q[0];
mcx q[18],q[17],q[12],q[10],q[9],q[8],q[1];
mcx q[18],q[17],q[15],q[12],q[10],q[2];
mcx q[18],q[17],q[15],q[12],q[10],q[1];
mcx q[18],q[17],q[12],q[10],q[3];
mcx q[18],q[17],q[14],q[13],q[12],q[10],q[9],q[0];
x q[10];
mcx q[18],q[17],q[16],q[12],q[3];
mcx q[18],q[17],q[16],q[12],q[0];
mcx q[18],q[17],q[16],q[15],q[12],q[2];
mcx q[18],q[17],q[16],q[15],q[12],q[1];
mcx q[18],q[17],q[16],q[15],q[12],q[0];
x q[11];
mcx q[18],q[12],q[11],q[10],q[3];
mcx q[18],q[15],q[12],q[11],q[10],q[2];
mcx q[18],q[15],q[12],q[11],q[10],q[1];
x q[10];
mcx q[18],q[12],q[11],q[10],q[9],q[2];
mcx q[18],q[15],q[14],q[12],q[11],q[10],q[1];
mcx q[18],q[14],q[12],q[11],q[10],q[9],q[7],q[0];
mcx q[18],q[15],q[12],q[11],q[10],q[8],q[7],q[0];
mcx q[18],q[15],q[13],q[12],q[11],q[10],q[8],q[0];
mcx q[18],q[16],q[12],q[11],q[3];
mcx q[18],q[16],q[12],q[11],q[2];
mcx q[18],q[16],q[12],q[11],q[1];
ccx q[17],q[11],q[5];
x q[13];
mcx q[15],q[14],q[13],q[9],q[0];
mcx q[15],q[13],q[9],q[8],q[0];
mcx q[18],q[17],q[13],q[12],q[10],q[9],q[8],q[0];
x q[10];
mcx q[18],q[13],q[12],q[11],q[10],q[9],q[8],q[0];
mcx q[18],q[16],q[15],q[13],q[12],q[11],q[8],q[0];
x q[16];
mcx q[18],q[16],q[13],q[12],q[11],q[9],q[8],q[0];
x q[14];
mcx q[16],q[14],q[10],q[9],q[1];
mcx q[15],q[14],q[9],q[1];
x q[10];
mcx q[18],q[17],q[15],q[14],q[13],q[12],q[10],q[0];
mcx q[18],q[17],q[15],q[14],q[12],q[10],q[0];
mcx q[18],q[17],q[14],q[12],q[10],q[9],q[7],q[0];
mcx q[18],q[17],q[14],q[12],q[10],q[9],q[1];
mcx q[18],q[14],q[12],q[11],q[10],q[9],q[1];
x q[10];
mcx q[16],q[14],q[13],q[10],q[9],q[0];
mcx q[16],q[15],q[14],q[13],q[10],q[0];
x q[16];
mcx q[18],q[17],q[16],q[14],q[13],q[12],q[9],q[0];
mcx q[18],q[17],q[16],q[15],q[14],q[13],q[12],q[0];
x q[7];
x q[10];
mcx q[18],q[12],q[11],q[10],q[9],q[8],q[7],q[0];
ccx q[13],q[7],q[1];
x q[13];
mcx q[15],q[14],q[9],q[7],q[0];
x q[8];
mcx q[18],q[16],q[15],q[12],q[11],q[8],q[1];
mcx q[18],q[16],q[12],q[11],q[9],q[8],q[1];
mcx q[18],q[17],q[16],q[12],q[9],q[8],q[1];
mcx q[18],q[17],q[16],q[12],q[9],q[8],q[0];
mcx q[18],q[17],q[16],q[15],q[12],q[8],q[1];
mcx q[18],q[17],q[16],q[15],q[12],q[8],q[0];
mcx q[18],q[17],q[15],q[13],q[12],q[10],q[8],q[0];
mcx q[18],q[17],q[15],q[12],q[10],q[8],q[1];
mcx q[15],q[9],q[8],q[7],q[0];
mcx q[15],q[9],q[8],q[0];
mcx q[18],q[16],q[15],q[12],q[11],q[8],q[7],q[0];
mcx q[18],q[17],q[16],q[15],q[12],q[8],q[7],q[0];
mcx q[18],q[17],q[16],q[12],q[9],q[8],q[7],q[0];
mcx q[18],q[17],q[12],q[10],q[9],q[8],q[7],q[0];
x q[17];
mcx q[18],q[16],q[12],q[11],q[9],q[8],q[7],q[0];
x q[11];
ccx q[14],q[8],q[2];
x q[15];
mcx q[17],q[15],q[11],q[10],q[2];
mcx q[17],q[15],q[11],q[10],q[7],q[0];
mcx q[17],q[16],q[15],q[14],q[11],q[7],q[0];
x q[7];
mcx q[17],q[15],q[11],q[10],q[8],q[7],q[0];
mcx q[17],q[15],q[14],q[11],q[10],q[0];
mcx q[17],q[15],q[14],q[13],q[11],q[10],q[0];
mcx q[17],q[16],q[15],q[11],q[2];
x q[14];
mcx q[17],q[15],q[14],q[11],q[10],q[1];
mcx q[17],q[16],q[15],q[14],q[11],q[1];
mcx q[17],q[16],q[15],q[11],q[8],q[7],q[0];
mcx q[17],q[16],q[15],q[13],q[11],q[8],q[0];
x q[16];
mcx q[17],q[15],q[13],q[11],q[10],q[8],q[0];
x q[10];
mcx q[16],q[15],q[10],q[8],q[1];
x q[7];
mcx q[16],q[15],q[14],q[10],q[7],q[0];
mcx q[17],q[15],q[14],q[11],q[10],q[7],q[0];
mcx q[17],q[15],q[14],q[11],q[7],q[0];
x q[13];
mcx q[16],q[15],q[10],q[8],q[7],q[0];
mcx q[16],q[15],q[10],q[7],q[0];
mcx q[16],q[15],q[13],q[10],q[8],q[0];
x q[16];
mcx q[17],q[16],q[15],q[14],q[13],q[11],q[0];
mcx q[17],q[16],q[15],q[13],q[11],q[0];
x q[16];
x q[13];
x q[11];
mcx q[18],q[16],q[15],q[12],q[11],q[2];
mcx q[18],q[16],q[15],q[12],q[11],q[1];
mcx q[18],q[15],q[14],q[13],q[12],q[11],q[10],q[0];
mcx q[18],q[16],q[15],q[14],q[13],q[12],q[11],q[0];
mcx q[18],q[15],q[12],q[11],q[10],q[8],q[1];
mcx q[18],q[15],q[12],q[11],q[8],q[1];
x q[16];
mcx q[18],q[16],q[15],q[14],q[12],q[11],q[1];
mcx q[18],q[16],q[15],q[14],q[12],q[11],q[0];
mcx q[18],q[16],q[15],q[14],q[12],q[11],q[7],q[0];
mcx q[18],q[16],q[15],q[12],q[11],q[7],q[0];
x q[11];
x q[17];
mcx q[18],q[17],q[16],q[15],q[14],q[12],q[1];
mcx q[18],q[17],q[16],q[15],q[13],q[12],q[8],q[0];
mcx q[18],q[17],q[16],q[15],q[14],q[12],q[7],q[0];
x q[7];
mcx q[18],q[17],q[16],q[15],q[12],q[7],q[0];
x q[17];
x q[9];
mcx q[17],q[16],q[11],q[9],q[2];
mcx q[17],q[16],q[11],q[9],q[8],q[7],q[0];
mcx q[17],q[16],q[11],q[9],q[7],q[0];
mcx q[17],q[16],q[14],q[11],q[9],q[7],q[0];
mcx q[17],q[16],q[14],q[11],q[9],q[1];
mcx q[17],q[16],q[13],q[11],q[9],q[8],q[0];
x q[16];
mcx q[16],q[13],q[10],q[9],q[8],q[0];
mcx q[16],q[10],q[9],q[8],q[1];
mcx q[16],q[10],q[9],q[8],q[0];
x q[8];
ccx q[15],q[9],q[3];
x q[7];
mcx q[16],q[10],q[9],q[8],q[7],q[0];
mcx q[16],q[14],q[10],q[9],q[7],q[0];
x q[10];
mcx q[17],q[11],q[10],q[9],q[2];
mcx q[17],q[14],q[11],q[10],q[9],q[1];
mcx q[17],q[14],q[11],q[10],q[9],q[0];
mcx q[17],q[11],q[10],q[9],q[8],q[7],q[0];
mcx q[17],q[14],q[11],q[10],q[9],q[7],q[0];
mcx q[17],q[11],q[10],q[9],q[8],q[0];
x q[8];
mcx q[17],q[13],q[11],q[10],q[9],q[8],q[0];
x q[14];
mcx q[17],q[16],q[14],q[13],q[11],q[9],q[0];
mcx q[17],q[14],q[13],q[11],q[10],q[9],q[0];
mcx q[17],q[14],q[13],q[11],q[9],q[0];
x q[14];
x q[11];
x q[16];
mcx q[18],q[16],q[14],q[12],q[11],q[9],q[1];
mcx q[18],q[16],q[14],q[12],q[11],q[9],q[0];
mcx q[18],q[16],q[14],q[13],q[12],q[11],q[9],q[0];
mcx q[18],q[14],q[13],q[12],q[11],q[10],q[9],q[0];
mcx q[18],q[12],q[11],q[10],q[9],q[8],q[1];
x q[17];
mcx q[18],q[17],q[16],q[14],q[12],q[9],q[1];
mcx q[18],q[17],q[16],q[13],q[12],q[9],q[8],q[0];
x q[14];
mcx q[18],q[16],q[14],q[12],q[11],q[9],q[7],q[0];
mcx q[18],q[17],q[16],q[14],q[12],q[9],q[7],q[0];
x q[12];
h q[12];
h q[6];
cx q[6], q[12];
h q[12];
h q[6];
ccx q[18],q[12],q[5];
ccx q[18],q[12],q[4];
ccx q[18],q[12],q[3];
ccx q[18],q[12],q[2];
ccx q[18],q[12],q[1];
x q[18];
cx q[18], q[12];
cx q[12], q[6];
cx q[18], q[12];
cx q[12], q[6];
h q[18];
h q[5];
cx q[5], q[18];
h q[18];
h q[5];
x q[16];
ccx q[16], q[18], q[4];
x q[16];
ccx q[16], q[18], q[4];
h q[18];
h q[3];
cx q[3], q[18];
h q[18];
h q[3];
cx q[18], q[10];
cx q[10], q[2];
cx q[18], q[10];
cx q[10], q[2];
h q[18];
h q[1];
cx q[1], q[18];
h q[18];
h q[1];