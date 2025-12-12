#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QTextEdit>
#include <QPixmap>
#include <QImage>
#include <QPainter>
#include <QFileDialog>
#include <QProgressBar>
#include <QMessageBox>
#include <QScrollArea>
#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QFile>
#include <QDataStream>
#include <QVector>
#include <QStatusBar>
#include <cmath>
#include <map>
#include <limits>

const double PI = 3.14159265358979323846;
const double NAN_VALUE = std::numeric_limits<double>::quiet_NaN();

// GSHHS Polygon structure
struct GSHHSPolygon
{
    int id;
    int n_points;
    int level;
    double west, east, south, north;
    QVector<double> lons;
    QVector<double> lats;
};

// Read 32-bit big-endian integer
int32_t readBigEndian32(QFile &file)
{
    unsigned char buffer[4];
    if (file.read((char *)buffer, 4) != 4)
        return 0;
    return ((int32_t)buffer[0] << 24) | ((int32_t)buffer[1] << 16) |
           ((int32_t)buffer[2] << 8) | (int32_t)buffer[3];
}

// Grid parameters structure
struct GridParams
{
    double lambda_0;
    double phi_0;
    double azimuth_sampling;
    double elevation_sampling;
};

// Get grid parameters for given SSD
GridParams getGridParams(double ssd)
{
    if (std::abs(ssd - 0.5) < 1e-6)
    {
        return {8.9142405037 * PI / 180, -8.9142405037 * PI / 180,
                0.000800524494 * PI / 180, 0.000800524494 * PI / 180};
    }
    else if (std::abs(ssd - 1.0) < 1e-6)
    {
        return {8.9138402398 * PI / 180, -8.9138402398 * PI / 180,
                0.001601048988 * PI / 180, 0.001601048988 * PI / 180};
    }
    else if (std::abs(ssd - 2.0) < 1e-6)
    {
        return {8.9130397083 * PI / 180, -8.9130397083 * PI / 180,
                0.003202097973 * PI / 180, 0.003202097973 * PI / 180};
    }
    throw std::invalid_argument("SSD must be 0.5, 1.0, or 2.0 km");
}

// FCI lat/lon to grid result
struct LatLonToGridResult
{
    QVector<double> rows;
    QVector<double> columns;
};

// FCI grid to lat/lon result
struct GridToLatLonResult
{
    QVector<double> lons;
    QVector<double> lats;
    QVector<bool> earth_mask;
};

long count_invisible;
long count_visible;


// Convert latitude/longitude to FCI grid row/column
LatLonToGridResult fci_latlon_to_grid(const QVector<double> &lats,
                                      const QVector<double> &lons,
                                      double ssd = 1.0)
{
    GridParams params = getGridParams(ssd);
    double r_eq = 6378.137;
    double r_pol = r_eq * (1 - 1.0 / 298.257223563);
    double h = 35786.4 + r_eq;

    QVector<double> result_rows, result_cols;

    for (int i = 0; i < lats.size(); ++i)
    {
        double lat_rad = lats[i] * PI / 180.0;
        double lon_rad = lons[i] * PI / 180.0;

        double cos_lat = std::cos(lat_rad);
        double sin_lat = std::sin(lat_rad);
        double cos_lon = std::cos(lon_rad);
        double sin_lon = std::sin(lon_rad);

        double c = 1.0 / std::sqrt(cos_lat * cos_lat +
                                   (r_pol / r_eq) * (r_pol / r_eq) * sin_lat * sin_lat);

        double x = r_eq * c * cos_lat * cos_lon;
        double y = r_eq * c * cos_lat * sin_lon;
        double z = r_pol * c * (r_pol / r_eq) * sin_lat;

        double x_sat = h;
        double dx = x - x_sat;
        double dy = y;
        double dz = z;

        double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
        double cos_angle = -dx / distance;
        double horizon_angle = std::asin(r_eq / h);
        //bool visible = cos_angle < std::cos(horizon_angle);
        double g = cos_lat * cos_lon;
        bool visible = g >= (1.0/distance);
        if(!visible)
            count_invisible += 1;
        else {
            count_visible += 1;
        }

        double rho_c = std::sqrt(dx * dx + dy * dy);
        double phi_s = std::atan(dz / rho_c);
        double lambda_s = std::atan(dy / dx);

        double col = (params.lambda_0 - lambda_s) / params.azimuth_sampling + 1.0;
        double row = (phi_s - params.phi_0) / params.elevation_sampling + 1.0;

        result_rows.push_back(visible ? row : NAN_VALUE);
        result_cols.push_back(visible ? col : NAN_VALUE);
    }

    return {result_rows, result_cols};
}

// Convert FCI grid row/column to latitude/longitude
GridToLatLonResult fci_grid_to_latlon(const QVector<double> &rows,
                                      const QVector<double> &cols,
                                      double ssd = 1.0)
{
    GridParams params = getGridParams(ssd);
    double r_eq = 6378.137;
    double r_pol = r_eq * (1 - 1.0 / 298.257223563);
    double h = 35786.4 + r_eq;
    double lambda_D = 0.0;

    QVector<double> result_lons, result_lats;
    QVector<bool> result_earth_mask;

    for (int i = 0; i < rows.size(); ++i)
    {
        double lambda_s = params.lambda_0 - (cols[i] - 1.0) * params.azimuth_sampling;
        double phi_s = params.phi_0 + (rows[i] - 1.0) * params.elevation_sampling;

        double cos_lambda_s = std::cos(lambda_s);
        double sin_lambda_s = std::sin(lambda_s);
        double cos_phi_s = std::cos(phi_s);
        double sin_phi_s = std::sin(phi_s);

        double s5 = h * h - r_eq * r_eq;
        double s4 = (r_eq / r_pol) * (r_eq / r_pol);

        double s_d_2 = (h * cos_lambda_s * cos_phi_s) * (h * cos_lambda_s * cos_phi_s) -
                       s5 * (cos_phi_s * cos_phi_s + s4 * sin_phi_s * sin_phi_s);

        bool is_earth = s_d_2 >= 0;

        if (!is_earth)
        {
            result_lons.push_back(NAN_VALUE);
            result_lats.push_back(NAN_VALUE);
            result_earth_mask.push_back(false);
            continue;
        }

        double s_d = std::sqrt(s_d_2);
        double s_n = (h * cos_lambda_s * cos_phi_s - s_d) /
                     (cos_phi_s * cos_phi_s + s4 * sin_phi_s * sin_phi_s);

        double s1 = h - s_n * cos_lambda_s * cos_phi_s;
        double s2 = -s_n * sin_lambda_s * cos_phi_s;
        double s3 = s_n * sin_phi_s;

        double s_xy = std::sqrt(s1 * s1 + s2 * s2);

        double lat = std::atan((s3 / s_xy) * s4);
        double lon = std::atan(s2 / s1) + lambda_D * PI / 180.0;

        result_lons.push_back(lon * 180.0 / PI);
        result_lats.push_back(lat * 180.0 / PI);
        result_earth_mask.push_back(true);
    }

    return {result_lons, result_lats, result_earth_mask};
}

// Read GSHHS binary file
QVector<GSHHSPolygon> readGSHHS(const QString &filepath)
{
    QVector<GSHHSPolygon> polygons;
    QFile file(filepath);

    if (!file.open(QIODevice::ReadOnly))
    {
        qWarning() << "Cannot open GSHHS file:" << filepath;
        return polygons;
    }

    while (!file.atEnd())
    {
        GSHHSPolygon poly;

        poly.id = readBigEndian32(file);
        poly.n_points = readBigEndian32(file);
        int flag = readBigEndian32(file);

        int west_i = readBigEndian32(file);
        int east_i = readBigEndian32(file);
        int south_i = readBigEndian32(file);
        int north_i = readBigEndian32(file);

        readBigEndian32(file); // area
        readBigEndian32(file); // area_full
        readBigEndian32(file); // container
        readBigEndian32(file); // ancestor

        poly.west = west_i * 1e-6;
        poly.east = east_i * 1e-6;
        poly.south = south_i * 1e-6;
        poly.north = north_i * 1e-6;

        poly.level = flag & 255; //(flag >> 24) & 255;

        for (int i = 0; i < poly.n_points; ++i)
        {
            int lon_i = readBigEndian32(file);
            int lat_i = readBigEndian32(file);

            double lon = lon_i * 1e-6;
            double lat = lat_i * 1e-6;

            // Normalize longitude from 0-360 to -180-180 range
            if (lon > 180.0)
            {
                lon -= 360.0;
            }

            poly.lons.push_back(lon);
            poly.lats.push_back(lat);
        }

        polygons.push_back(poly);
    }

    file.close();
    return polygons;
}

// Main window class
class ReferenceGridApp : public QMainWindow
{
    Q_OBJECT

public:
    ReferenceGridApp(QWidget *parent = nullptr) : QMainWindow(parent)
    {
        setWindowTitle("FCI Reference Grid & GSHHS Viewer");
        setGeometry(100, 100, 1200, 800);

        // Create central widget
        QWidget *centralWidget = new QWidget(this);
        QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);

        // Create tab widget
        QTabWidget *tabs = new QTabWidget(this);

        // Tab 1: FCI Grid Calculations
        QWidget *fciTab = createFCITab();
        tabs->addTab(fciTab, "FCI Grid");

        // Tab 2: GSHHS Coastlines
        QWidget *gshhsTab = createGSHHSTab();
        tabs->addTab(gshhsTab, "GSHHS Coastlines");

        mainLayout->addWidget(tabs);
        setCentralWidget(centralWidget);
    }

private slots:
    void onConvertCoordinates()
    {
        QVector<double> lats = {50.85, 51.51, 48.86};
        QVector<double> lons = {4.35, -0.13, 2.35};
        QStringList cities = {"Brussels", "London", "Paris"};

        QString result = "City to Grid Conversion (SSD=1.0):\n\n";

        auto grid_result = fci_latlon_to_grid(lats, lons, 1.0);

        for (int i = 0; i < cities.size(); ++i)
        {
            result += QString("%1 (%.2f°, %.2f°):\n")
                          .arg(cities[i])
                          .arg(lats[i])
                          .arg(lons[i]);
            result += QString("  Row: %.2f, Column: %.2f\n\n")
                          .arg(grid_result.rows[i])
                          .arg(grid_result.columns[i]);
        }

        outputTextEdit->setText(result);
    }

    void onLoadGSHHS()
    {
        QString filepath = "/home/hugo/EUMETCastTools/EUMETCastView/bin/gshhs2_3_7/gshhs_l.b";

        statusBar()->showMessage("Loading GSHHS file...");
        QApplication::processEvents();

        auto polygons = readGSHHS(filepath);

        if (polygons.empty())
        {
            QMessageBox::warning(this, "Error", "Failed to load GSHHS file");
            statusBar()->showMessage("Ready");
            return;
        }

        // Calculate statistics
        std::map<int, int> level_counts;
        int total_points = 0;

        for (const auto &poly : polygons)
        {
            level_counts[poly.level]++;
            total_points += poly.n_points;
        }

        // Display in table
        gshhsTable->setRowCount(0);
        gshhsTable->setColumnCount(3);
        gshhsTable->setHorizontalHeaderLabels({"Level", "Type", "Count"});

        QStringList level_names = {"", "Land", "Lake", "Island", "Pond"};
        int row = 0;

        for (auto &[level, count] : level_counts)
        {
            gshhsTable->insertRow(row);
            gshhsTable->setItem(row, 0, new QTableWidgetItem(QString::number(level)));
            gshhsTable->setItem(row, 1, new QTableWidgetItem(level_names[level]));
            gshhsTable->setItem(row, 2, new QTableWidgetItem(QString::number(count)));
            row++;
        }

        QString info = QString("GSHHS Statistics:\n\n"
                               "Total Polygons: %1\n"
                               "Total Points: %2\n\n"
                               "Level Distribution:\n")
                           .arg(polygons.size())
                           .arg(total_points);

        for (auto &[level, count] : level_counts)
        {
            if (level >= 1 && level <= 4)
            {
                info += QString("  Level %1 (%2): %3\n")
                            .arg(level)
                            .arg(level_names[level])
                            .arg(count);
            }
        }

        gshhsInfoLabel->setText(info);
        statusBar()->showMessage(QString("Loaded %1 polygons").arg(polygons.size()));
    }

    void onRenderMap()
    {
        QString filepath = "/home/hugo/EUMETCastTools/EUMETCastView/bin/gshhs2_3_7/gshhs_l.b";

        statusBar()->showMessage("Rendering map...");
        QApplication::processEvents();

        auto polygons = readGSHHS(filepath);

        if (polygons.empty())
        {
            QMessageBox::warning(this, "Error", "Failed to load GSHHS file");
            statusBar()->showMessage("Ready");
            return;
        }

        // Create image
        int width = 1024, height = 512;
        QImage image(width, height, QImage::Format_RGB32);
        image.fill(Qt::white);

        double map_west = -180.0, map_east = 180.0;
        double map_south = -90.0, map_north = 90.0;

        auto pixel_from_latlon = [&](double lat, double lon) -> QPair<int, int>
        {
            int x = (int)((lon - map_west) / (map_east - map_west) * width);
            int y = (int)((map_north - lat) / (map_north - map_south) * height);
            return {x, y};
        };

        auto is_valid_pixel = [&](int x, int y) -> bool
        {
            return x >= 0 && x < width && y >= 0 && y < height;
        };

        // Draw coastlines
        QPainter painter(&image);
        painter.setRenderHint(QPainter::Antialiasing);

        for (const auto &poly : polygons)
        {
            //if(poly.level > 4 || poly.level == 1)
            //    continue;
            QColor color;
            if (poly.level == 1)
                color = QColor(255, 0, 0);
            else if( poly.level == 2)
                color = QColor(0, 255, 0);
            else if (poly.level == 3)
                color = QColor(0, 0, 255);
            else  if (poly.level == 4)
                color = QColor(0, 255, 255);
            else    
                color = QColor(0, 0, 0);
                

            painter.setPen(color);

            for (int i = 1; i < poly.lons.size(); ++i)
            {
                double lon1 = poly.lons[i - 1];
                double lon2 = poly.lons[i];
                double lat1 = poly.lats[i - 1];
                double lat2 = poly.lats[i];

                // Normalize longitudes to handle wrap-around
                // If difference is > 180, adjust one coordinate
                double lon_diff = lon2 - lon1;
                if (lon_diff > 180.0)
                {
                    lon2 -= 360.0;
                }
                else if (lon_diff < -180.0)
                {
                    lon2 += 360.0;
                }

                auto p1 = pixel_from_latlon(lat1, lon1);
                auto p2 = pixel_from_latlon(lat2, lon2);

                // Only draw if both endpoints are within or near the image bounds
                if (is_valid_pixel(p1.first, p1.second) || is_valid_pixel(p2.first, p2.second))
                {
                    painter.drawLine(p1.first, p1.second, p2.first, p2.second);
                }
            }
        }
        painter.end();

        // Display image
        mapLabel->setPixmap(QPixmap::fromImage(image));
        image.save("gshhs_coastlines.png");
        statusBar()->showMessage("Map rendered and saved as gshhs_coastlines.png");
    }

    void onRenderPerspectiveView()
    {
        statusBar()->showMessage("Rendering perspective view...");
        QApplication::processEvents();

        const int SIZE = 11136;
        const double SSD = 1.0;
        count_invisible = 0;
        count_visible = 0;
        QImage image(SIZE, SIZE, QImage::Format_RGB32);
        image.fill(Qt::white);

        // Load GSHHS polygons
        QString filepath = "/home/hugo/EUMETCastTools/EUMETCastView/bin/gshhs2_3_7/gshhs_l.b";
        auto polygons = readGSHHS(filepath);

        if (polygons.empty())
        {
            QMessageBox::warning(this, "Error", "Failed to load GSHHS file");
            statusBar()->showMessage("Ready");
            return;
        }

        QPainter painter(&image);
        painter.setRenderHint(QPainter::Antialiasing);

        // For each polygon, convert all lat/lon points to grid coordinates
        for (const auto &poly : polygons)
        {
            // Convert lat/lon to grid coordinates
            auto result = fci_latlon_to_grid(poly.lats, poly.lons, SSD);

            //if(poly.level > 4 || poly.level == 1)
            //    continue;
            QColor color;
            if (poly.level == 1)
                color = QColor(255, 0, 0);
            else if( poly.level == 2)
                color = QColor(0, 255, 0);
            else if (poly.level == 3)
                color = QColor(0, 0, 255);
            else  if (poly.level == 4)
                color = QColor(0, 255, 255);
            else    
                color = QColor(0, 0, 0);


            painter.setPen(color);

            // Draw line segments between consecutive points
            for (int i = 1; i < result.rows.size(); ++i)
            {
                double row1 = result.rows[i - 1];
                double col1 = result.columns[i - 1];
                double row2 = result.rows[i];
                double col2 = result.columns[i];

                // Skip invalid points (NaN or out of bounds)
                if (std::isnan(row1) || std::isnan(col1) ||
                    std::isnan(row2) || std::isnan(col2))
                {
                    continue;
                }

                // Clamp to image bounds
                int x1 = (int)std::max(0.0, std::min((double)SIZE - 1, col1));
                int y1 = (int)std::max(0.0, std::min((double)SIZE - 1, row1));
                int x2 = (int)std::max(0.0, std::min((double)SIZE - 1, col2));
                int y2 = (int)std::max(0.0, std::min((double)SIZE - 1, row2));

                // Check for large jumps (polygon wrapping)
                if (std::abs(x2 - x1) < SIZE / 2 && std::abs(y2 - y1) < SIZE / 2)
                {
                    painter.drawLine(x1, y1, x2, y2);
                }
            }
        }
        painter.end();

        qDebug() << "Total invisible points:" << count_invisible;
        qDebug() << "Total visible points:" << count_visible;
        image.flip(Qt::Vertical);
        // Display and save image
        mapLabel->setPixmap(QPixmap::fromImage(image).scaledToWidth(512, Qt::SmoothTransformation));
        image.save("fci_perspective_view.png");
        statusBar()->showMessage("Perspective view rendered and saved as fci_perspective_view.png");
    }

private:
    QWidget *createFCITab()
    {
        QWidget *widget = new QWidget();
        QVBoxLayout *layout = new QVBoxLayout(widget);

        QLabel *titleLabel = new QLabel("FCI Level 1c Reference Grid Converter");
        QFont font = titleLabel->font();
        font.setPointSize(14);
        font.setBold(true);
        titleLabel->setFont(font);

        layout->addWidget(titleLabel);
        layout->addSpacing(20);

        QPushButton *convertBtn = new QPushButton("Convert Cities to Grid");
        connect(convertBtn, &QPushButton::clicked, this, &ReferenceGridApp::onConvertCoordinates);

        outputTextEdit = new QTextEdit();
        outputTextEdit->setReadOnly(true);
        outputTextEdit->setFont(QFont("Monospace"));

        layout->addWidget(convertBtn);
        layout->addWidget(outputTextEdit);

        return widget;
    }

    QWidget *createGSHHSTab()
    {
        QWidget *widget = new QWidget();
        QVBoxLayout *layout = new QVBoxLayout(widget);

        QLabel *titleLabel = new QLabel("GSHHS Coastline Data");
        QFont font = titleLabel->font();
        font.setPointSize(14);
        font.setBold(true);
        titleLabel->setFont(font);

        layout->addWidget(titleLabel);
        layout->addSpacing(10);

        // Buttons
        QHBoxLayout *buttonLayout = new QHBoxLayout();
        QPushButton *loadBtn = new QPushButton("Load GSHHS Data");
        QPushButton *renderBtn = new QPushButton("Render Map");
        QPushButton *perspectiveBtn = new QPushButton("Render Perspective View");
        connect(loadBtn, &QPushButton::clicked, this, &ReferenceGridApp::onLoadGSHHS);
        connect(renderBtn, &QPushButton::clicked, this, &ReferenceGridApp::onRenderMap);
        connect(perspectiveBtn, &QPushButton::clicked, this, &ReferenceGridApp::onRenderPerspectiveView);

        buttonLayout->addWidget(loadBtn);
        buttonLayout->addWidget(renderBtn);
        buttonLayout->addWidget(perspectiveBtn);
        buttonLayout->addStretch();

        layout->addLayout(buttonLayout);

        // Info and table
        QHBoxLayout *contentLayout = new QHBoxLayout();

        gshhsInfoLabel = new QLabel();
        gshhsInfoLabel->setAlignment(Qt::AlignTop);
        gshhsInfoLabel->setStyleSheet("QLabel { background-color: #f0f0f0; padding: 10px; }");

        gshhsTable = new QTableWidget();
        gshhsTable->setMaximumWidth(300);

        contentLayout->addWidget(gshhsInfoLabel);
        contentLayout->addWidget(gshhsTable);

        layout->addLayout(contentLayout);

        // Map display
        mapLabel = new QLabel();
        mapLabel->setAlignment(Qt::AlignCenter);
        mapLabel->setStyleSheet("QLabel { background-color: #e0e0e0; }");
        mapLabel->setMinimumHeight(400);

        layout->addWidget(new QLabel("Rendered Map:"));
        layout->addWidget(mapLabel);

        return widget;
    }

    QTextEdit *outputTextEdit;
    QLabel *gshhsInfoLabel;
    QTableWidget *gshhsTable;
    QLabel *mapLabel;
};

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    ReferenceGridApp window;
    window.show();

    return app.exec();
}

#include "reference_grid_qt1.moc"
