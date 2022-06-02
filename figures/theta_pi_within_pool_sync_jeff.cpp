/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2022 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include "genesis/genesis.hpp"
#include "genesis/population/formats/sam_variant_input_iterator.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <unordered_map>

using namespace genesis;
using namespace genesis::population;
using namespace genesis::utils;

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_INFO << "started";

    using VariantWindow = Window<genesis::population::Variant>;

    // ${GRENEDALF} diversity
    // --sync-path all_sims/${fname}
    // --window-width ${num_sites}
    // --pool-sizes ${pool_size}
    // --min-coverage 3
    // --measure all
    // --out-dir "pi"
    // --file-prefix "${fname%.*}-" > "pi_log/${fname%.*}.log"
    if( argc != 5 ) {
        throw std::runtime_error( "need input" );
    }

    auto const infile = std::string( argv[1] );
    auto const window_width = stoi( std::string( argv[2] ));
    auto const pool_size = stoi( std::string( argv[3] ));
    auto const prefix = std::string( argv[4] );

    auto it = make_variant_input_iterator_from_sync_file( infile );
    auto win_it = make_default_sliding_interval_window_iterator( it.begin(), it.end() );
    win_it.width( window_width );
    // win_it.stride( window_stride );

    // Tally up!
    double total_pi = 0;
    size_t total_snp = 0;
    size_t total_win = 0;
    LOG_TIME << "Starting iterator";
    for( auto const& window : win_it ) {

        auto sample_range = make_transform_range(
            []( VariantWindow::Entry const& entry ) -> BaseCounts const& {
                return entry.data.samples[0];
            },
            window.begin(), window.end()
        );

        total_pi += theta_pi_within_pool( sample_range.begin(), sample_range.end(), pool_size );
        total_snp += window.entry_count();
        ++total_win;
    }
    LOG_TIME << "Finished iterator";

    if( total_win != 1 ) {
        throw std::runtime_error(
            "window size should be way bigger than all contant of the simulations. "
            "something is wrong here."
        );
    }

    // 10_10_0.0_0_const_1_pop_theta_0.001-diversity.csv
//     chrom	start	end	1.snp_count	1.coverage_fraction	1.theta_pi_abs	1.theta_pi_rel	1.theta_watterson_abs	1.theta_watterson_rel	1.tajimas_d
// chr1	1	10000000	31954	0.004	19911.3851	0.445504655	19790.7264	0.442804994	0.0134326183
    std::ofstream csv_out( "pi/" + prefix + "diversity.csv" );
    csv_out << "chrom\tstart\tend\t1.snp_count\t1.coverage_fraction\t1.pi_within_abs\t1.pi_within_rel\n";
    csv_out << "chr1\t1\t" << window_width << "\t" << total_snp << "\t";
    csv_out << ( static_cast<double>(total_snp) / static_cast<double>(window_width) ) << "\t";
    csv_out << total_pi << "\t";
    csv_out << std::setprecision (15) << ( static_cast<double>(total_pi) / static_cast<double>(window_width) ) << "\n";

    return 0;
}
