

namespace mort
{
    static const int CELL_NBR = 3;

    double img_dis2(const double* v1, const double* v2, const numvec& box);

    numvec img_frac(const double* frac, const numvec& box);
  
    void get_frac(int natm3, const double* pos, const numvec& box, double* frac);

    struct cell_t
    {
    public:

        cell_t() {};
        
        void add_atom( int atomid, const numvec& imgcrd)
        {
            atomids.push_back(atomid);
            imgcrds.push_back(imgcrd[0]);
            imgcrds.push_back(imgcrd[1]);
            imgcrds.push_back(imgcrd[2]);
        }

        void add_nbrcell( int cellid )
        {
            nbrcells.push_back(cellid);
        }

        int size() const
        {
            return atomids.size();
        }

        int nnbr() const
        {
            return nbrcells.size();
        }
       
        vector<int> atomids;
        vector<int> nbrcells;
        vector<double> imgcrds;
    };

    class grid_t : public vector<cell_t>
    {
    public:

        grid_t() : m_box(3) {}

        void init( const numvec& boxsize, double cut )
        {
            m_box = boxsize;
       
            double cell_size = cut/CELL_NBR;
            m_dimsiz[0] = int( ceil(boxsize[0]/cell_size) + 0.01 );
            m_dimsiz[1] = int( ceil(boxsize[1]/cell_size) + 0.01 );
            m_dimsiz[2] = int( ceil(boxsize[2]/cell_size) + 0.01 );

            int ncell = m_dimsiz[0] * m_dimsiz[1] * m_dimsiz[2];
            resize(ncell);

            double cellen_x = m_box[0]/m_dimsiz[0];
            double cellen_y = m_box[1]/m_dimsiz[1];
            double cellen_z = m_box[2]/m_dimsiz[2];

            m_cut = max( makevec(cellen_x, cellen_y, cellen_z) )*CELL_NBR;

            m_offset[0] = 1;
            m_offset[1] = m_dimsiz[0];
            m_offset[2] = m_dimsiz[0] * m_dimsiz[1];
            
            m_offset_dot   = m_offset[0];
            m_offset_line  = m_offset[1];
            m_offset_layer = m_offset[2];

            for( int icell=0; icell < ncell; ++icell )
            {
                list_nbr_cell( icell );
            }
        }

        void list_nbr_cell( int icell )
        {
            int center = icell;
            
            int line_center = center;
            int line_id = line_center/m_offset_line;
           
            for( int inbr=0; inbr <= CELL_NBR; ++inbr )
            {
                int new_cell = line_center+inbr;

                if( new_cell/m_offset_line > line_id )
                {
                    new_cell -= m_offset_line;
                }

                assert( new_cell >= 0 );
                at( center ).add_nbrcell( new_cell );
            }

            int layer_center = center;
            int layer_id = layer_center/m_offset_layer;
            
            for( int jnbr=1; jnbr <= CELL_NBR; ++jnbr )
            {
                line_center = layer_center + jnbr*m_offset_line;

                if( line_center/m_offset_layer > layer_id )
                {
                    line_center -= m_offset_layer;
                }

                line_id = line_center/m_offset_line;
                for( int inbr= -CELL_NBR; inbr <= CELL_NBR; ++inbr )
                {
                    int new_cell = line_center+inbr;
                    
                    if( new_cell < 0 )
                    {
                        new_cell += m_offset_line;
                    }

                    if( new_cell/m_offset_line > line_id )
                    {
                        new_cell -= m_offset_line;
                    }
                    else if( new_cell/m_offset_line < line_id )
                    {
                        new_cell += m_offset_line;
                    }

                    assert( new_cell >= 0 );
                    at( center ).add_nbrcell( new_cell );
                }
            }

            for( int layer=1; layer <= CELL_NBR; ++layer )
            {
                int layer_center = center + layer*m_offset_layer;
 
                if( layer_center >= (int)size() )
                {
                    layer_center -= size();
                }
         
                int layer_id = layer_center/m_offset_layer;

                for( int line = -CELL_NBR; line <= CELL_NBR; ++line )
                {
                    line_center = layer_center + line*m_offset_line;

                    if( line_center < 0 )
                    {
                        line_center += m_offset_layer;
                    }

                    if( line_center/m_offset_layer < layer_id )
                    {
                        line_center += m_offset_layer;
                    }
                    else if( line_center/m_offset_layer > layer_id )
                    {
                        line_center -= m_offset_layer;
                    }

                    assert( line_center >= 0 );
                    line_id = line_center/m_offset_line;

                    for( int inbr = -CELL_NBR; inbr <= CELL_NBR; ++inbr )
                    {
                        int new_cell = line_center+inbr;

                        if( new_cell < 0 )
                        {
                            new_cell += m_offset_line;
                        }

                        if( new_cell/m_offset_line > line_id )
                        {
                            new_cell -= m_offset_line;
                        }
                        else if( new_cell/m_offset_line < line_id )
                        {
                            new_cell += m_offset_line;
                        }
                        
                        assert( new_cell >= 0 );
                        at( center ).add_nbrcell( new_cell );
                    }
                }
            }
        }

        void assign_atom(int natm3, double* frac)
        {
            int natom = natm3 / 3;
            for( int i=0; i < natom; ++i )
            {
                int i1 = int( (frac[3*i  ]+0.5) * m_dimsiz[0] );
                int i2 = int( (frac[3*i+1]+0.5) * m_dimsiz[1] );
                int i3 = int( (frac[3*i+2]+0.5) * m_dimsiz[2] );
                int cellid = i1*m_offset[0] + i2*m_offset[1] + i3*m_offset[2];
                at(cellid).add_atom( i, img_frac(&frac[3*i], m_box) );
            }
        }

        
    private:

        double m_cut;
        int m_dimsiz[3];
        int m_offset[3];
        int m_offset_dot;
        int m_offset_line;
        int m_offset_layer;
        numvec m_box;
    };

} // namespace mort

