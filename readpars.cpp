#include "readpars.h"

bool read_one_line(std::ifstream *ifs, std::string &name, std::string &val,
                   int &line_no)
{
    bool skip = true;
    char c;

    /* Skip all the comment lines */
    while (skip)
    {
        if (!ifs->good() || !ifs->get(c))
            return (false);

        if (c == '#') /* Comment line : skip it */
        {
            while (ifs->get(c) && (c != '\n'))
                ;
            ++line_no;
        }
        else
        {
            if ((c != '\n') && c) // skip empty lines
            {
                skip = false;
            }
        }
    }
    name += c;

    /* Read char until next space */
    while (ifs->good() && ifs->get(c) && (c != ' ') && (c != '\n') && (c != '\t'))
    {
        name += c;
    }

    if (!c || (c == '\n') || !ifs->good())
    {
        std::string error;
        std::stringstream ss;
        ss << "Error line " << line_no << ": Missing value for " << name;
        ss >> error;
        throw (std::invalid_argument(error));

    }

    /* Read separating spaces */
    while (ifs->good() && ifs->get(c) && ((c == ' ') || (c == '\t')))
        ;
    if (!c || (c == '\n') || !ifs->good())
    {
        std::string error;
        std::stringstream ss;
        ss << "Error line " << line_no << ": Missing value for " << name;
        ss >> error;
        throw (std::invalid_argument(error));

    }

    /* Read Value */
    if (c == '"') /* Read a string */
    {
        while (ifs->good() && ifs->get(c) && (c != '\n'))
        {
            if (c == '\\') /* Escape char */
            {
                if (!ifs->good() || !ifs->get(c) || (c == '\n'))
                {
                    /* Syntax error : \ at end of line */

                    std::string error;
                    std::stringstream ss;
                    ss << "Error line " << line_no << ": Escape character at end of line";
                    ss >> error;
                    throw (std::invalid_argument(error));

                }
                switch (c)
                {
                    case 'n':
                        c = '\n';
                        break;
                    case 't':
                        c = '\t';
                        break;
                }
                val += c;
            }
            else if (c == '"') /* End of string */
            {
                /* Read end of line */
                while (ifs->good() && ifs->get(c) && (c != '\n'))
                    ;
                ++line_no;
                return (true);
            }
            else
            {
                val += c;
            }
        }
        /* Syntax error : string not terminated */
        std::string error;
        std::stringstream ss;
        ss << "Error line " << line_no << ": string not terminated";
        ss >> error;
        throw (std::invalid_argument(error));

    }
    else /* Read Rest of the line until EOL or '#'*/
    {
        val += c;
        while (ifs->good() && ifs->get(c) && (c != '#') && (c != '\n'))
        {
            val += c;
        }
        if (!c || (c == '\n') || !ifs->good())
        {
            ++line_no;
            return (true);
        }
        /* Read end of line */
        while (ifs->good() && ifs->get(c) && (c != '\n'))
            ;
        ++line_no;
        return (true);
    }
    return (false);
}

/***********************************************/
/* Read the parameter file fname               */
/* Each line must be a parameter listed in par */
/* and in the same order                       */
/* Lines starting with a # are skipped         */
/***********************************************/
void read_pars(const char *fname, struct Parameter *par)
{
    std::ifstream *ifs;
    std::string Name, Val;
    int line_no = 0;
    ifs = new std::ifstream(fname);
    if (!ifs->good())
    {
        std::string error;
        error = " ERROR: Can't open file ";
        error+=fname;
        throw (std::invalid_argument(error));

    }

    while (par->name && read_one_line(ifs, Name, Val, line_no))
    {
        if (std::string(par->name) != Name)
        {

            std::string error;
            std::stringstream ss;
            cout << " ERROR: line " << line_no << " in " << fname
                 << " : expect parameter " << par->name;
            ss >> error;
            cout << error << endl;
            throw (std::invalid_argument::exception());

        }
        switch (par->type)
        {
            case T_INT:
                sscanf(Val.c_str(), "%d", (int *)(par->p));
                break;
            case T_DOUBLE:
                sscanf(Val.c_str(), "%lg", (double *)(par->p));
                break;
            case T_STRING:
                strcpy((char *)(par->p), Val.c_str());
                break;
            default:
                std::string error;
                std::stringstream ss;
                ss << " ERROR: : invalid parameter type " << par->type << "\n";
                ss >> error;
                cout << error << endl;
                throw (std::invalid_argument(error));

        }

        ++par;
        Name = Val = "";
    }
}

/***********************************************************/
/* Return the index of string name in the string array ptr */
/* Return : index of name in string array                  */
/*          -1 : name not found                            */
/***********************************************************/
int string_index(const char **ptr, const char *name)
{
    int i = 0;

    while (*ptr)
    {
        if (!strcmp(*ptr++, name))
            return (i);
        ++i;
    }
    return (-1);
}
