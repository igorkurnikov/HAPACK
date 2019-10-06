#include <sstream>
#include <boost/algorithm/string.hpp>
#include <common.hpp>
#include <object.hpp>

namespace mort
{
    using std::for_each;
    namespace ambfrc
    {

        namespace amoeba
        {

            void read_atom(istream& is, molecule_t& atomff, molecule_t& poleff)
            {
                int poleid, atomid;
                string name;

                is >> poleid >> atomid >> name;
                while (name[name.length() - 1] != '"')
                {
                    is >> name;
                }

                int element;
                double mass;
                is >> element >> mass;

                string line;
                std::getline(is, line);

                if (atomid > atomff.natom())
                {
                    assert( atomid==atomff.natom()+1 );
                    atom_t atom = atom_t::create( atomff, "");
                    atom.set_i(ELEMENT, element);
                    atom.set_d(MASS, mass);
                }

                assert( poleid==poleff.natom()+1 );
                atom_t pole = atom_t::create( poleff, "");
                pole.set_i(TYPEID, atomid);
            }

            void read_vdw(istream& is, molecule_t& atomff)
            {
                string sid;
                double rstar, depth;
                is >> sid >> rstar >> depth;

                atom_t atom(atomff, atoi(sid.c_str())-1);
                atom.set_s(NAME, sid);
                atom.set_d(RSTAR, rstar);
                atom.set_d(DEPTH, depth);

                assert( atoi(sid.c_str())==atom.get_i(ID) );

                string buff = " ";
                std::getline(is, buff);

                if (empty(buff))
                {
                    atom.set_d(VBUFF, 1.0);
                }
                else
                {
                    atom.set_d(VBUFF, 1.0 - atof(buff.c_str()));
                }

            }

            void read_bond(istream& is, molecule_t& mol)
            {
                int id1, id2;
                double force, equil;
                is >> id1 >> id2 >> force >> equil;

                atom_t atm1(mol, id1-1);
                atom_t atm2(mol, id2-1);
                bond_t bond = bond_t::create(atm1, atm2);

                bond.set_d(FORCE, force);
                bond.set_d(EQUIL, equil);
            }

            void read_angl(istream& is, molecule_t& atomff)
            {
                int id1, id2, id3;
                double force, equil;
                is >> id1 >> id2 >> id3 >> force >> equil;
                atom_t atm1(atomff, id1-1);
                atom_t atm2(atomff, id2-1);
                atom_t atm3(atomff, id3-1);

                // assert( has_bond( atm1, atm2 ) && has_bond( atm2, atm3 ) );
                angl_t angl = angl_t::create(atm1, atm2, atm3);
                angl.set_d(FORCE, force);
                angl.set_d(EQUIL, equil);

                string line;
                std::getline(is, line);

                if (!empty(line))
                {
                    numvec addition = makevec(0.0, 0.0);
                    std::istringstream ls(line);
                    ls >> addition[0] >> addition[1];
                    angl.set_v(ADDITION, addition);
                }
            }

            void read_urey(istream& is, molecule_t& atomff)
            {
                if (!atomff.angl_begin()->has_v(UREY))
                {
                    for_each(atomff.angl_begin(), atomff.angl_end(),
                            vparm_setter1(UREY, makevec(0.0, 0.0)));
                }

                int id1, id2, id3;
                numvec urey(2);
                is >> id1 >> id2 >> id3 >> urey[0] >> urey[1];

                atom_t atm1(atomff, id1 - 1);
                atom_t atm2(atomff, id2 - 1);
                atom_t atm3(atomff, id3 - 1);
                angl_t angl = angl_t::get(atm1, atm2, atm3);

                angl.set_v(UREY, urey);
            }

            void read_strbnd(istream& is, molecule_t& atomff)
            {
                int id;
                numvec strbnd(3);
                is >> id >> strbnd[0] >> strbnd[1] >> strbnd[2];

                atom_t atom(atomff, id-1);

                if (!atom.has_v(STRBND))
                {
                    for_each(atomff.atom_begin(), atomff.atom_end(),
                            vparm_setter1(STRBND, makevec(0.0, 0.0, 0.0)));
                }

                atom.set_v(STRBND, strbnd);
            }

            void read_oops(istream& is, molecule_t& atomff)
            {
                int id1, id2;
                double force;
                is >> id1 >> id2 >> force;

                atom_t atm1(atomff, id1-1);
                atom_t atm2(atomff, id2-1);

                assert( bond_t::has(atm1, atm2) );
                bond_t bond = bond_t::get(atm1, atm2);

                if (!bond.has_v(OPBEND))
                {
                    for_each(atomff.bond_begin(), atomff.bond_end(),
                            vparm_setter1(OPBEND, makevec(0.0, 0.0)));
                }

                numvec opbend = bond.get_v(OPBEND);
                if (atom_1st(bond) == atom_2nd(bond))
                {
                    opbend[0] = force;
                    opbend[1] = force;
                }
                else if (atom_1st(bond) == atm1)
                {
                    opbend[0] = force;
                }
                else
                {
                    opbend[1] = force;
                }

                bond.set_v(OPBEND, opbend);
                assert( bond.get_v(OPBEND,opbend) );
            }

            void read_tors(istream& is, molecule_t& mol)
            {
                atomvec_t atoms;
                for (int i = 0; i < 4; ++i)
                {
                    int id;
                    is >> id;
                    atoms.push_back(atom_t(mol, id-1));
                }

                for (int i = 0; i < 3; ++i)
                {
                    int period;
                    double force, phase;
                    is >> force >> phase >> period;

                    dihe_t tors = dihe_t::create(atoms[0],atoms[1],atoms[2],atoms[3],period);
                    tors.set_d(FORCE, force);
                    tors.set_d(EQUIL, phase * M_PI / 180.0);
                }
            }

            void read_pitors(istream& is, molecule_t& mol)
            {
                int id1, id2;
                double force;
                is >> id1 >> id2 >> force;

                atom_t atm1(mol, id1 - 1);
                atom_t atm2(mol, id2 - 1);
                bond_t bond = bond_t::get(atm1, atm2);
                bond.set_d(PITORS, force);
            }

            void read_tortor(istream& is, molecule_t& mol)
            {
                atomvec_t atoms;
                for (int i = 0; i < 5; i++)
                {
                    int id;
                    is >> id;
                    atoms.push_back(atom_t(mol, id - 1));

                }

                shared_ptr<vector<double> > pfuncs(new vector<double> ());

                int slice_1, slice_2;
                is >> slice_1 >> slice_2;
                for (int j = 0; j < slice_1; j++)
                {
                    for (int k = 0; k < slice_2; k++)
                    {
                        double phi1, phi2, func, dfd1, dfd2, dfda;
                        is >> phi1;
                        is >> phi2;
                        is >> func;
                        is >> dfd1;
                        is >> dfd2;
                        is >> dfda;

                        pfuncs->push_back(phi1);
                        pfuncs->push_back(phi2);
                        pfuncs->push_back(func);
                        pfuncs->push_back(dfd1);
                        pfuncs->push_back(dfd2);
                        pfuncs->push_back(dfda);
                    }
                }

                morf_t tor2 = create_tor2(atoms);
                tor2.set_a(TOR2FUNC, pfuncs);
            }

            void read_mpole(istream& is, molecule_t& poleff)
            {
                std::string line;
                std::getline(is, line);
				 
                vector<string> items;
                boost::split(items, line, is_any_of(" "), token_compress_on);

                numvec pole(9);
                is >> pole[0] >> pole[1] >> pole[2];
                is >> pole[3] >> pole[4] >> pole[5];
                is >> pole[6] >> pole[7] >> pole[8];

                if( !( items.size() == 5 || items.size() == 6) ) throw std::runtime_error(" invalid multipole line " + line );  

                if (items.size() == 5)
                {
                    int id1 = std::abs(atoi(items[1].c_str()));
                    int id2 = std::abs(atoi(items[2].c_str()));
                    int id3 = std::abs(atoi(items[3].c_str()));

                    if (id2 > 0 && id3 > 0)
                    {
                        atom_t atm1(poleff, id1 - 1);
                        atom_t atm2(poleff, id2 - 1);
                        atom_t atm3(poleff, id3 - 1);
                        angl_t angl = angl_t::create(atm2, atm1, atm3);

                        if (items[2].find('-') != string::npos
                                || items[3].find('-') != string::npos)
                        {
                            angl.set_i(TYPEID, (int) BISECTOR);
                        }
                        else
                        {
                            angl.set_i(TYPEID, (int) MONO);
                        }

                        angl.set_d(PCHG, atof(items[4].c_str()));
                        angl.set_v(POLE, pole);
                    }
                }
                else
                {
                    atom_t a1(poleff, std::abs(atoi(items[1].c_str())) - 1);
                    atom_t a2(poleff, std::abs(atoi(items[2].c_str())) - 1);
                    atom_t a3(poleff, std::abs(atoi(items[3].c_str())) - 1);
                    atom_t a4(poleff, std::abs(atoi(items[4].c_str())) - 1);

                    impr_t oops = impr_t::create(a2,a3,a1,a4,2);
                    oops.set_d(PCHG, atof(items[5].c_str()));
                    oops.set_v(POLE, pole);
                }
            }

            void read_polar(istream& is, molecule_t& poleff)
            {
                string sid;
                double polar;
                is >> sid >> polar;

                int iid = atoi(sid.c_str());
                assert( iid <= poleff.natom() );

                atom_t atom(poleff, iid - 1);
                atom.set_s(NAME, sid);
                atom.set_d(POLAR, polar);

                string line;
                std::getline(is, line);

                int id2;
                std::istringstream ls(line);
                while (ls >> id2)
                {
                    atom_t next(poleff, id2 - 1);
                    bond_t::frcget(atom, next);
                }
            }

        } // namespace amoeba

    } // namespace ambfrc

    using namespace ambfrc;

    void read_amoeba_frc(istream& is, molecule_t& atomff, molecule_t& poleff)
    {
//		std::cerr << " read_amoeba_frc pt 1 " << std::endl;

        string keyw = " ";

        atomff.set_s(NAME, "amoeba-atom");
        poleff.set_s(NAME, "amoeba-pole");

        while (is >> keyw)
        {
//            std::cout << " read_amoeba_frc(): reading keyw: " << keyw << std::endl;

            if (keyw[0] == '#')
            {
                is.ignore(MAX_LINE_WIDTH, '\n');
            }
            else if (keyw == "atom")
            {
                //std::cout << "read amoeba atom" << std::endl;
                amoeba::read_atom(is, atomff, poleff);
            }
            else if (keyw == "vdw")
            {
                //std::cout << "read amoeba vdw" << std::endl;
                amoeba::read_vdw(is, atomff);
            }
            else if (keyw == "bond")
            {
                //std::cout << "read amoeba bond" << std::endl;
                amoeba::read_bond(is, atomff);
            }
            else if (keyw == "angle")
            {
                //std::cout << "read amoeba angl" << std::endl;
                amoeba::read_angl(is, atomff);
            }
            else if (keyw == "torsion")
            {
                //std::cout << "read amoeba torsion" << std::endl;
                amoeba::read_tors(is, atomff);
            }
            else if (keyw == "strbnd")
            {
                //std::cout << "read amoeba strbnd" << std::endl;
                amoeba::read_strbnd(is, atomff);
            }
            else if (keyw == "opbend")
            {
                //std::cout << "read amoeba opbend" << std::endl;
                amoeba::read_oops(is, atomff);
            }
            else if (keyw == "ureybrad")
            {
                //std::cout << "read amoeba ureybrad" << std::endl;
                amoeba::read_urey(is, atomff);
            }
            else if (keyw == "pitors")
            {
                //std::cout << "read amoeba pitors" << std::endl;
                amoeba::read_pitors(is, atomff);
            }
            else if (keyw == "tortors")
            {
                //std::cout << "read amoeba tortor" << std::endl;
                amoeba::read_tortor(is, atomff);
            }
            else if (keyw == "multipole")
            {
                //std::cout << "read amoeba multipole" << std::endl;
                amoeba::read_mpole(is, poleff);
            }
            else if (keyw == "polarize")
            {
                //std::cout << "read amoeba polar" << std::endl;
                amoeba::read_polar(is, poleff);
            }
            else if (keyw == "biotype")
            {
                string line;
                std::getline(is, line);
            }
            else
            {
                string value = next_word(is);
                is.ignore(MAX_LINE_WIDTH, '\n');
                atomff.set_d(keyw.c_str(), atof(value.c_str()));
            }
        }
    }

} // namespace mort


