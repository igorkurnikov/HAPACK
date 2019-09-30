template <class T1, class T2>
T2
static_cast_with_check(T1 __ptr)
{
  g_assert(__ptr);
  return static_cast<T2>(__ptr);
}


// 
// macros
// 
#define G_TYPE_CHECK_INSTANCE_TYPE_FROM_UNQUOTE_NAME(ptr, tname)        \
    G_TYPE_CHECK_INSTANCE_TYPE(ptr, g_type_from_name(#tname))

// ui definition helper macros
#define GENERIC_STATIC_SIGNAL0_HANDLER(static_function_name,		\
				       member_function_name,		\
				       class_name,			\
				       signal_class)			\
    static void								\
    static_function_name(signal_class* __signal_class__ptr,		\
                         gpointer __user_data_ptr)			\
    {									\
        g_assert(__user_data_ptr);                                      \
        static_cast_with_check<gpointer, class_name*>(__user_data_ptr)-> \
            member_function_name();                                     \
    }								

// generic gtk widget event handler
#define GENERIC_STATIC_EVENT1_HANDLER(static_function_name,		\
				      member_function_name,		\
				      class_name,			\
				      event_class)			\
    static gboolean							\
    static_function_name(GtkWidget* __widget__ptr,			\
                         event_class* __event_ptr,			\
                         gpointer __user_data_ptr)			\
    {									\
        g_assert(__user_data_ptr && __event_ptr);                       \
        return static_cast_with_check<gpointer, class_name*>(__user_data_ptr)-> \
            member_function_name(__event_ptr);			   	\
    }								

// MEMBER function part
#define DECLARE_CLASS_ACTION_LISTENER(action_name)	\
    void on_##action_name##_action_listener();

#define DECLARE_CLASS_MENU_HINT_LISTENER(menu_name, hint_type)  \
    void on_##menu_name##_menu_##hint_type##_listener();

#define DECLARE_CLASS_MENU_HINTS_LISTENER(menu)		\
    DECLARE_CLASS_MENU_HINT_LISTENER(menu, select)	\
    DECLARE_CLASS_MENU_HINT_LISTENER(menu, deselect)

// static part
#define STATIC_ACTION_LISTENER_NAME(action_name)	\
    S_on_##action_name##_action_listener

#define STATIC_MENU_HINT_LISTENER_NAME(menu_name, hint)	\
    S_on_##menu_name##_##hint##_listener

#define STATIC_MENU_DESELECT_LISTENER_NAME(menu_name)	\
    S_on_##menu_name##_deselect_listener


#define DEFINE_STATIC_ACTION_LISTENER(class_name, action_name)		\
    GENERIC_STATIC_SIGNAL0_HANDLER(STATIC_ACTION_LISTENER_NAME(action_name), \
                                   on_##action_name##_action_listener,	\
                                   class_name,				\
                                   GtkAction) 

#define DEFINE_STATIC_MENU_HINT_LISTENER(class_name, menu_name, hint)	\
    GENERIC_STATIC_SIGNAL0_HANDLER(STATIC_MENU_HINT_LISTENER_NAME(menu_name, \
                                                                  hint), \
                                   on_##menu_name##_menu_##hint##_listener, \
                                   class_name,				\
                                   GtkItem) 

#define DEFINE_STATIC_MENU_HINTS_LISTENER(class, menu)          \
    DEFINE_STATIC_MENU_HINT_LISTENER(class, menu, select)	\
    DEFINE_STATIC_MENU_HINT_LISTENER(class, menu, deselect)	

#define CONNECT_ACTION_HINT_SIGNAL(action_path_ptr, callback)		\
    {									\
        g_assert(_M_ui_manager_ptr);					\
        GtkWidget* menu_item_ptr =					\
            gtk_ui_manager_get_widget(_M_ui_manager_ptr,                \
                                      action_path_ptr);			\
        g_assert(menu_item_ptr);                                        \
        g_signal_connect(menu_item_ptr,					\
                         "select",					\
                         G_CALLBACK(STATIC_MENU_HINT_LISTENER_NAME(callback, \
                                                                   select)), \
                         this);                                         \
        g_signal_connect(menu_item_ptr,					\
                         "deselect",					\
                         G_CALLBACK(STATIC_MENU_HINT_LISTENER_NAME(callback, \
                                                                   deselect)), \
                         this);                                         \
    }

#define DEFINE_MENU_HINT_FUNCTION(name, select_hint)    \
    void                                                \
    main_window::on_##name##_menu_select_listener()     \
    {                                                   \
        m_status.push(select_hint);                     \
    }                                                   \
    void                                                \
    main_window::on_##name##_menu_deselect_listener()   \
    {                                                   \
        m_status.pop();                                 \
    }

#define DEFINE_MENU_ITEM(menu, hint)				\
    DEFINE_STATIC_ACTION_LISTENER(main_window, menu)		\
    DEFINE_STATIC_MENU_HINTS_LISTENER(main_window, menu)        \
    DEFINE_MENU_HINT_FUNCTION(menu, hint)

#define SHOW_UI_WIDGET(path)					\
    {								\
        g_assert(m_pmanager);					\
        GtkWidget* item_ptr =                                   \
            gtk_ui_manager_get_widget(m_pmanager, path);        \
        g_assert(item_ptr);                                     \
        gtk_widget_show(item_ptr);                              \
    }

#define HIDE_UI_WIDGET(path)					\
    {								\
        g_assert(_M_ui_manager_ptr);                            \
        GtkWidget* item_ptr =                                   \
            gtk_ui_manager_get_widget(_M_ui_manager_ptr, path); \
        g_assert(item_ptr);                                     \
        gtk_widget_hide(item_ptr);                              \
    }
