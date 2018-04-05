#include "matcl-sqlite-cpp/builder/Select.h"
#include "matcl-sqlite-cpp/builder/update.h"
#include "matcl-sqlite-cpp/builder/insert.h"
#include "matcl-sqlite-cpp/builder/delete.h"
#include <iostream>

using namespace matcl::sql;

int syntax_test(int argc, char** argv) // syntax testing
{
    (void)argc;
    (void)argv;

    column id("id");
    column name("name");
    column number("number");
    table t1("t1");
    table t2("t2");
    table t3("t3", "t3_alias");
    table t4("t4");
    id.set_table(t2);
    name.set_table(t3);
    join_source js = t1.join(t2, number == id).join(t3.join(t4, id < 6, join_type::left), column_list()(id)(number));
    expression e = id > 3 || 50 <= id << 2 && concat(name, std::string("abc")) || is_not(name, null);
    std::string s = select().distinct()(id)(name).from( js ).where(e).to_str();
    
    std::cout << s << std::endl;

    s = select().where(id > 3 && is(bind_param("abc"), null) || call("foo", wildcard::all)).order_by(id)
        .to_str();

    std::cout << s << std::endl;

    s = select().order_by(ord_expression_list()(ord_expression(number,ordering::desc))
                        (ord_expression(collate(id,collation::rtrim)))).to_str();

    std::cout << s << std::endl;

    s = select().from(t3).order_by(number).limit(10, 20).to_str();
    
    std::cout << s << std::endl;

    s = select()(name).group_by(id, call("count", number) > 13).to_str();
    
    std::cout << s << std::endl;

    s = select()(id, "cze")(t4).intersect(select().from(t4)).limit(666).to_str();
    
    std::cout << s << std::endl;
    
    s = insert().columns( column_list()(id)(name) ).into(t2).values( expression_list()(1)("abc") ).to_str();

    std::cout << s << std::endl;

    s = insert().into( t3 ).to_str();

    std::cout << s << std::endl;

    s = update().on_conflict(conflict_behaviour::replace).set(assignment_list()(id = 3)(name = "bob")).tab(t1).to_str();
    
    std::cout << s << std::endl;

    s = delete_stmt().tab(t4).where(id == 4).to_str();
    
    std::cout << s << std::endl;

    s = update().on_conflict(conflict_behaviour::ignore).tab(table("main.matrix_table"))
        .set(column("mat_string") = bind_param("update_value"))
        .where(column("name") == bind_param("update_mat_name")).to_str();

    std::cout << s << std::endl;

    std::cin >> s;
    return 0;
}

