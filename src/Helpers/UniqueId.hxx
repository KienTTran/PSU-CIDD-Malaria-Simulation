/*
 * UniqueId.hxx
 * 
 * Define a thread-safe singleton pattern class that allows for unique identifiers
 * to be returned for the lifecycle of a single simulation. This is useful for 
 * some of the deeper model identification and validation to be done.
 */
#ifndef UNIQUEID_HXX
#define UNIQUEID_HXX

class UniqueId {
    private:
    long current_id;

        // Constructor
        UniqueId() { current_id = 0; }

        // Deconstructor
        ~UniqueId() = default;


    public:
        // Not supported by singleton
        UniqueId(UniqueId const&) = delete;

        // Not suported by singleton
        void operator=(UniqueId const&) = delete;

        // Get a reference to the object
        static UniqueId& get_instance() {
            static UniqueId instance;
            return instance;
        }

        // Return a unique id value
        long get_uid() {
            return ++current_id;
        }
};

#endif