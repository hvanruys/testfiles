cmake_minimum_required(VERSION 3.16)

project(ReferenceGridQt VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable automatic handling of Qt features
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Find Qt6 packages (change to Qt5 if needed)
find_package(Qt6 REQUIRED COMPONENTS
    Core
    Gui
    Widgets
)

# Alternative for Qt5 (uncomment if using Qt5 instead of Qt6):
# find_package(Qt5 REQUIRED COMPONENTS
#     Core
#     Gui
#     Widgets
# )

# Add the executable
add_executable(reference_grid_qt1
    reference_grid_qt1.cpp
)

# Link Qt libraries
target_link_libraries(reference_grid_qt1 PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
)

# Alternative for Qt5 (uncomment if using Qt5):
# target_link_libraries(reference_grid_qt1 PRIVATE
#     Qt5::Core
#     Qt5::Gui
#     Qt5::Widgets
# )

# Platform-specific settings
if(WIN32)
    set_target_properties(reference_grid_qt1 PROPERTIES
        WIN32_EXECUTABLE TRUE
    )
endif()

if(APPLE)
    set_target_properties(reference_grid_qt1 PROPERTIES
        MACOSX_BUNDLE TRUE
    )
endif()

if(UNIX AND NOT APPLE)
    # Linux-specific settings
    set_target_properties(reference_grid_qt1 PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    
    # Optional: Add desktop entry and icon support
    # Uncomment and customize if you want to install desktop integration
    # install(FILES reference_grid_qt1.desktop
    #     DESTINATION share/applications
    # )
    # install(FILES reference_grid_qt1.png
    #     DESTINATION share/icons/hicolor/256x256/apps
    # )
endif()

# Install target
install(TARGETS reference_grid_qt1
    BUNDLE DESTINATION .
    RUNTIME DESTINATION bin
)