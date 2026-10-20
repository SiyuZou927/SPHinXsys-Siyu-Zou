#ifndef ENTRY_BATCH_CONFIG_H
#define ENTRY_BATCH_CONFIG_H
#include <stdexcept>
#include <sstream>

struct EntryBatchConfig
{
    std::string mode, wave_type;
    Real speed, theta_deg, center_x, center_y, restart_time, duration, sample_dt, vtp_dt;
    int restart_step;

    EntryBatchConfig()
    {
        XmlEngine xml("entry_configuration", "parameters");
        xml.loadXmlFile("./input/entry_case.xml");
        auto element = xml.root_element_.element_begin("Entry");
        if (element == xml.root_element_.element_end())
            throw std::runtime_error("Missing Entry configuration element");
        xml.getRequiredAttributeValue(element, "mode", mode);
        xml.getRequiredAttributeValue(element, "wave_type", wave_type);
        xml.getRequiredAttributeValue(element, "speed", speed);
        xml.getRequiredAttributeValue(element, "theta_deg", theta_deg);
        xml.getRequiredAttributeValue(element, "center_x", center_x);
        xml.getRequiredAttributeValue(element, "center_y", center_y);
        xml.getRequiredAttributeValue(element, "restart_time", restart_time);
        xml.getRequiredAttributeValue(element, "restart_step", restart_step);
        xml.getRequiredAttributeValue(element, "duration", duration);
        xml.getRequiredAttributeValue(element, "sample_dt", sample_dt);
        xml.getRequiredAttributeValue(element, "vtp_dt", vtp_dt);
        if ((mode != "prepare" && mode != "new_entry") ||
            (wave_type != "regular" && wave_type != "focused") || speed <= 0 ||
            theta_deg <= 0 || theta_deg >= 90 || duration <= 0 || sample_dt <= 0 ||
            vtp_dt <= 0 || restart_step <= 0)
            throw std::runtime_error("Invalid entry configuration");
    }
};
#endif
