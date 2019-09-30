#ifndef DDUTILS_MORT_CREATE_H
#define DDUTILS_MORT_CREATE_H


namespace mort
{
    class molecule_t;

    bool create_ribbon( molecule_t& mol );

    void render_ribbon( const molecule_t& mol );

} // namespace mort

#endif
