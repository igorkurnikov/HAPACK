#ifndef MORT_ENERGY_CTRLPARM_HPP
#define MORT_ENERGY_CTRLPARM_HPP

#include <vector>
#include <string>

namespace mort
{
    using std::vector;
    using std::string;

    struct ctrlparm_t
    {
        ctrlparm_t(const string& str);

        void parse(const vector<string>& terms);

	int igb;

	double cut;

	double rgbmax;

	double offset;
    };

} // namespace mort


#endif
