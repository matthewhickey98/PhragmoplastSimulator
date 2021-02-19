#ifndef FNAME_H
#define FNAME_H

class FrameName {
public:
    FrameName() : index_(-1)
    {
        prefix_[0] = suffix_[0] = 0;
    };
    FrameName(const char *prefix, const char *suffix);

    char *next_name();
    char *current_name()
    {
        return (buffer_);
    }
    char *current_name_other_suffix(const char *suffix);
    void inc()
    {
        ++index_;
    }
    void set_index(int index)
    {
        index_ = index;
    }
    int index()
    {
        return (index_);
    };
    void set_prefix(const char *prefix);
    void set_suffix(const char *suffix);

protected:
    static const int MAX_FILE_SIZE = 1024;

    char prefix_[MAX_FILE_SIZE];
    char suffix_[MAX_FILE_SIZE];
    char buffer_[2 * MAX_FILE_SIZE];
    char buffer2_[2 * MAX_FILE_SIZE];
    int index_;
};

#endif
