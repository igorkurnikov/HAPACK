#ifndef MORT_OBJECT_ATOMEXPR_HPP
#define MORT_OBJECT_ATOMEXPR_HPP

#include <vector>
#include <boost/shared_ptr.hpp>


namespace mort
{
    using std::vector;
    using boost::shared_ptr;

    class morf_t;

    namespace atomexpr
    {
        enum env_e 
        { 
            HYBRID, 
            DEGREE, 
            VALENCE, 
            HYDROGEN, 
            NEBR
        };

        enum level_e   
        { 
            SINGLE, 
            CONN, 
            SURPRISE, 
            COMMA, 
            SEMICOLON 
        };

        enum bracket_e 
        { 
            OPEN, 
            CLOSE 
        };

        class node_i
        {
        public:

            node_i( int level );

            virtual ~node_i();

            virtual bool match( const morf_t& atom ) const = 0;

            int level();

            node_i* father();
            
            void set_father( node_i* father );            
            
        private:

            int m_level;

            node_i* m_father;
            
        };

        class iparm_node : public node_i
        {
          public:

            iparm_node(const hashid_t& parmid, int ivalue );

            virtual bool match(const morf_t& a) const;

          private:

            hashid_t m_parmid;
            int m_ivalue;
        };

        class env_node : public node_i
        {
          public:

            env_node(const hashid_t& env, int value );
            
            virtual ~env_node();
            
            virtual bool match( const morf_t& atom ) const;
            
          private:

            hashid_t m_env;
            
            int m_value;
        };    

        class logic_node : public node_i
        {
          public:

            logic_node( int level, const hashid_t& logic );
            
            ~logic_node();
            
            void push( const shared_ptr< node_i >& node );
            
            void pop();
            
            shared_ptr< node_i > top() const;

            int logic() const;

            virtual bool match( const morf_t& atom ) const;

          private:

            std::vector< shared_ptr< node_i > > m_children;

            hashid_t m_logic;   
        };


        class ring_node : public node_i
        {
          public:

            ring_node( int arom, int size );
            
            virtual ~ring_node();
         
            virtual bool match( const morf_t& atom  ) const;

          private:

            int m_size;
            
            int m_arom;
            
        };

        
        class true_node : public node_i
        {
        public:

            true_node() : node_i( SINGLE ) {}
            
            virtual ~true_node() {}

            virtual bool match( const morf_t& ) const { return true; }
        };

    } // namespace atomexpr

    class atomexpr_t
    {
      public:

        atomexpr_t();
            
        virtual ~atomexpr_t();
            
        void change_to(int level, const hashid_t& type);
            
        void insert(const shared_ptr< atomexpr::node_i >& node);
            
        bool match(const morf_t& atom) const;

      private:

        atomexpr_t(const atomexpr_t& rhs);
            
      private:

        shared_ptr< atomexpr::node_i > m_root;

        atomexpr::node_i* m_curt;
    };

} // namespace mort

#endif


