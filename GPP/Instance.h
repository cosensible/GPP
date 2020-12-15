#ifndef INSTANCE_H
#define INSTANCE_H

#include <string>

struct Instance {
	std::string insPath = "E://SMART_Work//GPP//Deploy//graphs//";
	std::string solPath = "E://SMART_Work//GPP//Deploy//solutions//";

	std::string insName = "fe_ocean";
	int partNum = 64;

	int nodeNum = 0;
	int edgeNum = 0;
	int isolated = 0;
	int maxDegree = 5000;
};

#endif // INSTANCE_H
